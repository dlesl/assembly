extern crate bio;
extern crate itertools;
#[macro_use]
extern crate failure;
#[cfg_attr(test, macro_use)]
extern crate gb_io;
#[macro_use]
extern crate log;

use bio::alignment::sparse;
use bio::alphabets::dna::revcomp;
use gb_io::seq::{Seq, Topology};

use itertools::Itertools;
use std::cmp;
use std::collections::{HashMap, HashSet};
use std::iter::Iterator;
use std::mem;
use std::ops::{Index, IndexMut};
use std::borrow::Borrow;
pub use crate::MatchIdx::*;

const MAX_PATH_LEN: usize = 100;
const MAX_PATHS: usize = 100_000;
const MAX_CPATHS: usize = 10_000;

#[derive(Debug, Fail)]
pub enum AssemblyError {
    #[fail(display = "Assembly too complex, MAX_PATHS exceeded")]
    Complexity,
}

#[derive(Hash, Eq, PartialEq, Ord, PartialOrd, Debug, Clone, Copy)]
pub enum MatchIdx {
    Idx(u32),
    IdxRc(u32),
}

impl MatchIdx {
    pub fn index(self) -> usize {
        match self {
            Idx(i) | IdxRc(i) => i as usize,
        }
    }
    pub fn flipped(self) -> MatchIdx {
        match self {
            Idx(i) => IdxRc(i),
            IdxRc(i) => Idx(i),
        }
    }
}

/// Describes a match between two fragments (a and b)
/// 0 - index into fragment a
/// 1 - index into fragment b
/// 2 - length of match
/// In the case of an `IdxRc`, the indices above refer to the reverse
/// complement of the sequence.
#[derive(Eq, PartialEq, Debug, Clone, Ord, PartialOrd, Hash)]
pub struct Match(pub u32, pub u32, pub u32);

/// A data structure allowing simple lookups of all matches between
/// a group of sequences. For every sequence with index `i`, the array
/// contains a HashMap for `Idx(i)` and `IdxRc(i)`, containing the
/// matches to the sequence and its reverse complement, respectively.
///
/// The HashMap's key is the index of the matching sequence, and the
/// value is a list of matches between the two (see `Match`).
///
/// `MatchArray` stores a lot of redundant information, however a
/// `Match` is very small and this simplifies its traversal when
/// searching for products later.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MatchArray(Vec<HashMap<MatchIdx, Vec<Match>>>);

impl MatchArray {
    pub fn new(num_entries: usize) -> MatchArray {
        MatchArray(vec![HashMap::new(); num_entries * 2]) // double length to make room for Rc's
    }
    /// Two-way insert - a -> b and b -> a
    pub fn insert(&mut self, from: MatchIdx, to: MatchIdx, matches: Vec<Match>) {
        let swap_matches = |matches: &[Match]| {
            matches
                .iter()
                .map(|&Match(pos_a, pos_b, length)| Match(pos_b, pos_a, length))
                .collect::<Vec<_>>()
        };
        self[to]
            .entry(from)
            .or_insert_with(Vec::new)
            .extend(swap_matches(&matches));
        self[from]
            .entry(to)
            .or_insert_with(Vec::new)
            .extend(matches);
    }
    pub fn keys(&self) -> impl Iterator<Item = MatchIdx> {
        self.keys_fwd().chain(self.keys_rc())
    }
    pub fn keys_fwd(&self) -> impl Iterator<Item = MatchIdx> {
        let n = self.0.len() / 2;
        (0..n).map(|i| Idx(i as u32))
    }
    pub fn keys_rc(&self) -> impl Iterator<Item = MatchIdx> {
        let n = self.0.len() / 2;
        (0..n).map(|i| IdxRc(i as u32))
    }
}

impl Index<MatchIdx> for MatchArray {
    type Output = HashMap<MatchIdx, Vec<Match>>;
    fn index(&self, rhs: MatchIdx) -> &HashMap<MatchIdx, Vec<Match>> {
        match rhs {
            Idx(i) => &self.0[i as usize],
            IdxRc(i) => &self.0[i as usize + (self.0.len() / 2)],
        }
    }
}

impl IndexMut<MatchIdx> for MatchArray {
    fn index_mut(&mut self, rhs: MatchIdx) -> &mut HashMap<MatchIdx, Vec<Match>> {
        match rhs {
            Idx(i) => &mut self.0[i as usize],
            IdxRc(i) => {
                let len = self.0.len();
                &mut self.0[i as usize + len / 2]
            }
        }
    }
}

pub trait SeqType /*: ToOwned*/ {
    fn len(&self) -> usize {
        self.seq().len()
    }
    fn seq(&self) -> &[u8];
}

impl<'a> SeqType for &'a [u8] {
    fn seq(&self) -> &[u8] {
        self
    }
}

impl SeqType for Seq {
    fn seq(&self) -> &[u8] {
        &self.seq
    }
}

impl<'a> SeqType for &'a Seq {
    fn seq(&self) -> &[u8] {
        &self.seq
    }
}

// inspired by: https://github.com/BjornFJohansson/pydna/blob/py3dev/pydna/assembly.py

pub fn find_homology<T>(seqs: &[T], limit: usize) -> MatchArray
where
    T: SeqType,
{
    let seqs: Vec<_> = seqs.iter().map(|s| s.seq().to_ascii_uppercase()).collect();
    let mut res = MatchArray::new(seqs.len());

    // hash the sequences
    let seqs: Vec<_> = seqs
        .iter()
        .map(|s| {
            let hashed = sparse::hash_kmers(s.as_slice(), limit);
            let rc = revcomp(s);
            (s, hashed, rc)
        })
        .collect();

    // check every sequence against its own rc
    for (i, (_, hashed, rc)) in seqs.iter().enumerate() {
        let matches = sparse::find_kmer_matches_seq1_hashed(hashed, rc, limit);
        let matches = consolidate_matches(matches, limit);
        if !matches.is_empty() {
            assert!(res[Idx(i as u32)].get(&IdxRc(i as u32)).is_none());
            assert!(res[IdxRc(i as u32)].get(&Idx(i as u32)).is_none());
            res.insert(Idx(i as u32), IdxRc(i as u32), matches);
        }
    }

    for ((ia, (a, a_hashed, _)), (ib, (b, _, b_rc))) in seqs.iter().enumerate().tuple_combinations()
    {
        // a against b (implies a_rc against b_rc)
        let matches = sparse::find_kmer_matches_seq1_hashed(a_hashed, b, limit);
        let matches = consolidate_matches(matches, limit);
        // we also store the relationships between all the sequences' reverse
        // complements - this closure 'flips' the matches.
        let reverse_matches = |matches: &[Match]| {
            matches
                .iter()
                .map(|&Match(pos_a, pos_b, length)| {
                    Match(
                        a.len() as u32 - pos_a - length,
                        b.len() as u32 - pos_b - length,
                        length,
                    )
                })
                .collect()
        };
        let matches_rc: Vec<_> = reverse_matches(&matches);
        assert!(res[Idx(ia as u32)].get(&Idx(ib as u32)).is_none());
        assert!(res[IdxRc(ia as u32)].get(&IdxRc(ib as u32)).is_none());
        if !matches.is_empty() {
            res.insert(Idx(ia as u32), Idx(ib as u32), matches);
            res.insert(IdxRc(ia as u32), IdxRc(ib as u32), matches_rc);
        }

        // a against b_rc (implies a_rc against b)
        let matches = sparse::find_kmer_matches_seq1_hashed(a_hashed, b_rc, limit);
        let matches = consolidate_matches(matches, limit);
        let matches_rc = reverse_matches(&matches);
        assert!(res[Idx(ia as u32)].get(&IdxRc(ib as u32)).is_none());
        assert!(res[IdxRc(ia as u32)].get(&Idx(ib as u32)).is_none());
        if !matches.is_empty() {
            res.insert(Idx(ia as u32), IdxRc(ib as u32), matches);
            res.insert(IdxRc(ia as u32), Idx(ib as u32), matches_rc);
        }
    }
    res
}

/// A node in a pathway through the "graph" of fragments
/// 0 - the index of fragment
/// 1 - describes how this fragment connects to the previous,
/// see `Match` for details
#[derive(Eq, PartialEq, Debug, Clone, Ord, PartialOrd, Hash)]
pub struct Node(pub MatchIdx, pub Match);

pub type Path = Vec<Node>;

#[derive(Eq, PartialEq, Debug, Clone)]
pub struct Products {
    pub linear: Vec<Path>,
    pub circular: Vec<Path>,
}

pub struct ProductsBuilder<'a, T: SeqType> {
    linear: HashSet<Path>,
    circular: HashSet<Path>,
    seqs: &'a [T],
}

impl<'a, T: SeqType> ProductsBuilder<'a, T> {
    fn new(seqs: &'a [T]) -> ProductsBuilder<'a, T> {
        ProductsBuilder {
            linear: HashSet::new(),
            circular: HashSet::new(),
            seqs,
        }
    }
    fn push_circular(&mut self, path: &Path) {
        if !self.circular.contains(&rotate(&invert(&path, self.seqs))) {
            self.circular.insert(rotate(&path));
        }
    }
    fn push_linear(&mut self, path: Path) {
        if !self.linear.contains(&invert(&path, self.seqs)) {
            self.linear.insert(path);
        }
    }
    fn products(self) -> Products {
        Products {
            linear: self.linear.into_iter().collect(),
            circular: self.circular.into_iter().collect(),
        }
    }
}

/// "rotates" a circular path so that it begins at the sequence with the lowest
/// index - this allows us to check for equality.
fn rotate(p: &Path) -> Path {
    assert!(!p.is_empty());
    let split_at = p.iter().enumerate().min_by_key(|&(_, x)| x).unwrap().0;
    p[split_at..]
        .iter()
        .chain(p[..split_at].iter())
        .cloned()
        .collect()
}

/// "inverts" a path, swapping every sequence for its reverse complement
fn invert<T: SeqType>(p: &Path, seqs: &[T]) -> Path {
    let mut res = Vec::with_capacity(p.len());
    for (Node(a_idx, a), Node(b_idx, _)) in p
        .iter()
        .rev()
        .cycle()
        .skip(p.len() - 1)
        .take(p.len() + 1)
        .tuple_windows()
    {
        res.push(Node(
            b_idx.flipped(),
            Match(
                seqs[a_idx.index()].len() as u32 - a.1 - a.2,
                seqs[b_idx.index()].len() as u32 - a.0 - a.2,
                a.2,
            ),
        ));
    }
    if (res[0].1).2 == 0 {
        // relinearise
        res[0].1 = Match(0, 0, 0);
    }
    res
}

/// Searches the `MatchArray` returned by `find_homology` for linear and circular products
pub fn find_products<T: SeqType>(
    matches: &MatchArray,
    seqs: &[T],
) -> Result<Products, AssemblyError> {
    let mut pb = ProductsBuilder::new(seqs);
    for k in matches.keys_fwd() {
        search(
            matches,
            seqs,
            &mut pb,
            &mut vec![Node(k, Match(0, 0, 0))], // Dummy match, to be replaced later when we "close" the path.
            0,
        )?;
    }
    Ok(pb.products())
}

/// A recursive depth-first search, we limit ourselves to a depth of `MAX_PATH_LEN` to avoid stack overflow
fn search<'a, T: SeqType>(
    matches: &MatchArray,
    seqs: &[T],
    pb: &mut ProductsBuilder<'a, T>,
    path: &mut Path,
    cur_pos: u32,
) -> Result<(), AssemblyError> {
    if path.len() > MAX_PATH_LEN || pb.linear.len() > MAX_PATHS || pb.circular.len() > MAX_CPATHS {
        return Err(AssemblyError::Complexity);
    }
    let start_at = path.last().unwrap().0;
    if path.len() > 1 {
        // We have a unique linear path
        pb.push_linear(path.clone());
        // Now check if we can close the path
        if let Some(matches) = matches[start_at].get(&path[0].0) {
            for m in matches {
                // check if 5' -> 3', and if it matches before the overlap to path[1]
                if m.0 > cur_pos && m.1 < (path[1].1).0 {
                    // Close the path (replace the dummy match we inserted at the beginning)
                    let mut path = path.clone();
                    path[0].1 = m.clone();
                    pb.push_circular(&path);
                }
            }
        }
    }
    for (idx, ms) in &matches[start_at] {
        if !path.iter().any(|&Node(ref n, _)| n == idx) {
            // each node can only be included once, maybe this excludes some possibilities?
            for m in ms {
                if m.0 > cur_pos &&
                    // check if this sequence has something to contribute (i.e. we're not matching its end)
                    seqs[idx.index()].len() as u32 - m.1 - m.2 > 0
                {
                    path.push(Node(*idx, m.clone()));
                    search(matches, seqs, pb, path, m.1)?;
                    path.pop();
                }
            }
        }
    }
    Ok(())
}

/// Returns a linear product if the first `Match` is zero length, otherwise circular
pub fn extract_product_bare<T: SeqType>(path: &Path, seqs: &[T]) -> Vec<u8> {
    assert!(path.len() > 1);
    let mut res = Vec::new();
    {
        let mut extend = |idx: &MatchIdx, a: u32, b: u32| match *idx {
            Idx(i) => {
                res.extend(&seqs[i as usize].seq()[a as usize..b as usize]);
            }
            IdxRc(i) => {
                let seq = seqs[i as usize].seq();
                res.extend(revcomp(
                    &seq[seq.len() - b as usize..seq.len() - a as usize],
                ));
            }
        };
        for (&Node(ref idx, Match(_, a, _)), &Node(_, Match(b, _, _))) in
            path.iter().tuple_windows()
        {
            extend(idx, a, b);
        }
        let circular = (path[0].1).2 != 0;
        if circular {
            let &Node(ref idx, Match(_, a, _)) = path.last().unwrap();
            let &Node(_, Match(b, _, _)) = path.first().unwrap();
            extend(idx, a, b);
        } else {
            let &Node(ref idx, Match(_, a, _)) = path.last().unwrap();
            extend(idx, a, seqs[idx.index()].len() as u32);
        }
    }
    res
}

pub fn product_len<T: SeqType>(path: &Path, seqs: &[T]) -> u32 {
    let mut res = 0;
    for (i, Node(idx, Match(_, to, _))) in path.iter().enumerate() {
        let next_from = path
            .get(i + 1)
            .map(|&Node(_, Match(f, _, _))| f)
            .unwrap_or_else(|| {
                if (path[0].1).2 == 0 {
                    // linear
                    seqs[idx.index()].len() as u32
                } else {
                    // circular
                    (path[0].1).0
                }
            });
        res += next_from - to;
    }
    res
}

/// Extracts products retaining annotations
pub fn extract_product_seq<T: Borrow<Seq>>(path: &Path, seqs: &[T]) -> Seq {
    assert!(path.len() > 1);
    let circular = (path[0].1).2 != 0;
    let mut res = Seq::empty();
    if circular {
        res.topology = Topology::Circular;
    }
    {
        let mut extend = |idx: &MatchIdx, a: u32, b: u32, len: u32| {
            let (mut seq, truncate) = match *idx {
                Idx(i) => {
                    let seq = seqs[i as usize].borrow();
                    let end = cmp::min(i64::from(b + len), seq.len());
                    let res = seq.extract_range(i64::from(a), end);
                    (res, end - i64::from(b))
                }
                IdxRc(i) => {
                    let seq = seqs[i as usize].borrow();
                    let b_rc = seq.len() - i64::from(b);
                    let end = cmp::max(b_rc, 0);
                    let mut res = seq.extract_range(end, seq.len() - i64::from(a)).revcomp();
                    let res_len = res.seq.len();
                    res.seq.truncate(res_len - (b_rc - end) as usize);
                    (res, b_rc - end)
                }
            };
            let features = mem::replace(&mut seq.features, Vec::new());
            for f in features {
                let offset = res.len();
                res.features.push(
                    seq.relocate_feature(f, offset)
                        .expect("Relocating feature failed"),
                );
            }
            res.seq
                .extend(&seq.seq[..seq.seq.len() - truncate as usize]);
        };
        for (&Node(ref idx, Match(_, a, _)), &Node(_, Match(b, _, len))) in
            path.iter().tuple_windows()
        {
            extend(idx, a, b, len);
        }
        if circular {
            let &Node(ref idx, Match(_, a, _)) = path.last().unwrap();
            let &Node(_, Match(b, _, len)) = path.first().unwrap();
            extend(idx, a, b, len);
        } else {
            let &Node(ref idx, Match(_, a, _)) = path.last().unwrap();
            extend(idx, a, seqs[idx.index()].borrow().len() as u32, 0);
        }
    }
    if circular {
        let features = mem::replace(&mut res.features, Vec::new());
        for mut f in features {
            match res.wrap_location(f.location) {
                Ok(l) => {
                    f.location = l;
                    res.features.push(f);
                }
                Err(e) => {
                    warn!("Error processing feature: {}", e);
                }
            }
        }
    }
    res
}

/// Consolidate adjacent matches to one long match, storing its length
fn consolidate_matches(matches: Vec<(u32, u32)>, k: usize) -> Vec<Match> {
    let k = k as u32;
    let mut res = Vec::with_capacity(matches.len());
    for (a, b) in matches {
        if res.is_empty() {
            res.push(Match(a, b, k));
        } else {
            let Match(last_a, last_b, count) = res.last().unwrap().clone();
            if last_a + count + 1 - k == a && last_b + count + 1 - k == b {
                match res.last_mut() {
                    Some(v) => {
                        v.2 += 1;
                    }
                    _ => unreachable!(),
                }
            } else {
                res.push(Match(a, b, k));
            }
        }
    }
    res
}

#[cfg(test)]
mod tests {
    use super::Match;
    use super::Node;
    use super::Path;
    use super::{Idx, IdxRc};
    use std::collections::HashSet;
    use std::iter::FromIterator;

    #[test]
    fn test1() {
        let test = vec![&b"gatc"[..], &b"agag"[..], &b"ttaa"[..]];
        let res = super::find_homology(&test, 4);
        assert_eq!(res.0.iter().filter(|r| !r.is_empty()).count(), 4);
        assert_eq!(res[Idx(0)].get(&IdxRc(0)), Some(&vec![Match(0, 0, 4)]));
        assert_eq!(res[Idx(2)].get(&IdxRc(2)), Some(&vec![Match(0, 0, 4)]));
        assert_eq!(res[IdxRc(0)].get(&Idx(0)), Some(&vec![Match(0, 0, 4)]));
        assert_eq!(res[IdxRc(2)].get(&Idx(2)), Some(&vec![Match(0, 0, 4)]));
    }

    #[test]
    fn consolidate_matches() {
        assert_eq!(
            super::consolidate_matches(vec![(0, 0)], 20),
            vec![Match(0, 0, 20)]
        );
        assert_eq!(
            super::consolidate_matches(vec![(0, 0), (1, 1)], 20),
            vec![Match(0, 0, 21)]
        );
        assert_eq!(
            super::consolidate_matches(vec![(0, 0), (1, 1), (2, 4)], 20),
            vec![Match(0, 0, 21), Match(2, 4, 20)]
        );
        assert_eq!(
            super::consolidate_matches(vec![(0, 0), (1, 1), (3, 2)], 20),
            vec![Match(0, 0, 21), Match(3, 2, 20)]
        );
        assert_eq!(
            super::consolidate_matches(vec![(0, 0), (1, 1), (3, 3)], 20),
            vec![Match(0, 0, 21), Match(3, 3, 20)]
        );
        assert_eq!(
            super::consolidate_matches(vec![(4, 6), (5, 7), (6, 10)], 20),
            vec![Match(4, 6, 21), Match(6, 10, 20)]
        );
        let test = &[&b"abcdefghijklm"[..], &b"defghi"[..]][..];
        let res = super::find_homology(&test, 4);
        assert_eq!(res.0.len(), 4);
        assert_eq!(res[Idx(0)].get(&Idx(1)), Some(&vec![Match(3, 0, 6)]));
        assert_eq!(res[Idx(1)].get(&Idx(0)), Some(&vec![Match(0, 3, 6)]));
        assert_eq!(res[IdxRc(0)].get(&IdxRc(1)), Some(&vec![Match(4, 0, 6)]));
        assert_eq!(res[IdxRc(1)].get(&IdxRc(0)), Some(&vec![Match(0, 4, 6)]));
    }

    #[test]
    fn test_assembly() {
        let one = &b"TCTGTATGAACGGTCTGGTCtttgccgaccgcacgccgcaTCCAGCGCTGACGGA"[..];
        let two = &b"TCCAGCGCTGACGGAgacgccgctgcgcgatcagttcacccgtgcaccgctgGATAACGACATTGGCGTA"[..];
        let three = &b"GATAACGACATTGGCGTAcaccgcatccggcgcggattggcctgaactgccagctggcgcaggtagcagagTCTGTATGAACGGTCTGGTC"[..];
        let seqs = &[one, two, three][..];
        let matches = super::find_homology(seqs, 15);
        let products = super::find_products(&matches, seqs).unwrap();
        assert_eq!(
            products.circular,
            vec![[
                Node(Idx(0), Match(71, 0, 20)),
                Node(Idx(1), Match(40, 0, 15)),
                Node(Idx(2), Match(52, 0, 18))
            ]]
        );
        // The order is random because of HashSet so compare like this:
        assert_eq!(
            products.linear.into_iter().collect::<HashSet<Path>>(),
            HashSet::from_iter(vec![
                vec![
                    Node(Idx(0), Match(0, 0, 0)),
                    Node(Idx(1), Match(40, 0, 15)),
                    Node(Idx(2), Match(52, 0, 18))
                ],
                vec![Node(Idx(0), Match(0, 0, 0)), Node(Idx(1), Match(40, 0, 15))],
                vec![Node(Idx(2), Match(0, 0, 0)), Node(Idx(0), Match(71, 0, 20))],
                vec![Node(Idx(1), Match(0, 0, 0)), Node(Idx(2), Match(52, 0, 18))],
                vec![
                    Node(Idx(2), Match(0, 0, 0)),
                    Node(Idx(0), Match(71, 0, 20)),
                    Node(Idx(1), Match(40, 0, 15))
                ],
                vec![
                    Node(Idx(1), Match(0, 0, 0)),
                    Node(Idx(2), Match(52, 0, 18)),
                    Node(Idx(0), Match(71, 0, 20))
                ]
            ])
        );
    }
    #[test]
    fn test_invert_linear() {
        let one = &b"TCTGTATGAACGGTCTGGTCtttgccgaccgcacgccgcaTCCAGCGCTGACGGA"[..];
        let two = &b"TCCAGCGCTGACGGAgacgccgctgcgcgatcagttcacccgtgcaccgctgGATAACGACATTGGCGTA"[..];
        let three = &b"GATAACGACATTGGCGTAcaccgcatccggcgcggattggcctgaactgccagctggcgcaggtagcagagTCTGTATGAACGGTCTGGTC"[..];
        let seqs = &[one, two, three][..];
        let input = vec![
            Node(IdxRc(0), Match(0, 0, 0)),
            Node(IdxRc(2), Match(35, 0, 20)),
            Node(IdxRc(1), Match(73, 0, 18)),
        ];
        let output = vec![
            Node(Idx(1), Match(0, 0, 0)),
            Node(Idx(2), Match(52, 0, 18)),
            Node(Idx(0), Match(71, 0, 20)),
        ];
        assert_eq!(super::invert(&input, seqs), output);
    }
    #[test]
    fn test_invert_circular() {
        use super::rotate;
        use super::MatchIdx::*;
        use super::Node;
        let one = &b"TCTGTATGAACGGTCTGGTCtttgccgaccgcacgccgcaTCCAGCGCTGACGGA"[..];
        let two = &b"TCCAGCGCTGACGGAgacgccgctgcgcgatcagttcacccgtgcaccgctgGATAACGACATTGGCGTA"[..];
        let three = &b"GATAACGACATTGGCGTAcaccgcatccggcgcggattggcctgaactgccagctggcgcaggtagcagagTCTGTATGAACGGTCTGGTC"[..];
        let seqs = &[one, two, three][..];
        let input = vec![
            Node(IdxRc(0), Match(55, 0, 15)),
            Node(IdxRc(2), Match(35, 0, 20)),
            Node(IdxRc(1), Match(73, 0, 18)),
        ];
        let output = vec![
            Node(Idx(0), Match(71, 0, 20)),
            Node(Idx(1), Match(40, 0, 15)),
            Node(Idx(2), Match(52, 0, 18)),
        ];
        assert_eq!(rotate(&super::invert(&input, seqs)), output);
    }

    #[test]
    fn test_extract_bare() {
        let one = b"abcdef";
        let two = b"defghi";
        let three = b"ghijklabc";
        let seqs = &[&one[..], &two[..], &three[..]][..];
        assert_eq!(
            super::extract_product_bare(
                &vec![
                    Node(Idx(0), Match(0, 0, 0)),
                    Node(Idx(1), Match(3, 0, 3)),
                    Node(Idx(2), Match(3, 0, 3))
                ],
                seqs
            ),
            b"abcdefghijklabc".to_vec()
        );
        assert_eq!(
            super::extract_product_bare(
                &vec![
                    Node(Idx(0), Match(6, 0, 3)),
                    Node(Idx(1), Match(3, 0, 3)),
                    Node(Idx(2), Match(3, 0, 3))
                ],
                seqs
            ),
            b"abcdefghijkl".to_vec()
        );
        let one = "attgcc".as_bytes();
        let two = "ttaggc".as_bytes();
        let three = "taannnatt".as_bytes();
        let seqs = &[one, two, three][..];
        assert_eq!(
            super::extract_product_bare(
                &vec![
                    Node(Idx(0), Match(0, 0, 0)),
                    Node(IdxRc(1), Match(3, 0, 3)),
                    Node(Idx(2), Match(3, 0, 3))
                ],
                seqs
            ),
            b"attgcctaannnatt".to_vec()
        );
    }
    #[test]
    fn test_extract_seq() {
        use super::*;
        use gb_io::seq::*;
        let make_seq = |positions: Vec<Location>| Seq {
            seq: "gatcgatgat".into(),
            topology: Topology::Linear,
            features: positions
                .into_iter()
                .map(|l| Feature {
                    location: l,
                    kind: feature_kind!(""),
                    qualifiers: Vec::new(),
                })
                .collect(),
            ..Seq::empty()
        };
        let seqs = vec![
            make_seq(vec![Location::simple_range(0, 10)]),
            make_seq(vec![Location::simple_range(0, 10)]),
        ];
        let check = |v1, v2| {
            let res = super::extract_product_seq(&v1, &seqs);
            println!("{}", String::from_utf8_lossy(&res.seq));
            assert_eq!(
                res.features.into_iter().map(|f| f.location).collect::<Vec<_>>(),
                v2
            );
        };
        check(
            vec![Node(Idx(0), Match(0, 0, 0)), Node(Idx(1), Match(7, 0, 3))],
            vec![Location::simple_range(0, 10), Location::simple_range(7, 17)],
        );
        check(
            vec![Node(Idx(0), Match(0, 0, 0)), Node(IdxRc(1), Match(7, 0, 3))],
            vec![
                Location::simple_range(0, 10),
                Location::Complement(Box::new(Location::simple_range(7, 17))),
            ],
        );
        check(
            vec![Node(Idx(0), Match(7, 0, 3)), Node(Idx(1), Match(7, 0, 3))],
            vec![
                Location::simple_range(0, 10),
                Location::Join(vec![
                    Location::simple_range(7, 14),
                    Location::simple_range(0, 3),
                ]),
            ],
        );
    }
}
