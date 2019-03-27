extern crate env_logger;
extern crate log;

#[macro_use]
extern crate clap;
extern crate assembly;
extern crate failure;
#[macro_use]
extern crate gb_io;

use std::fs;
use std::fs::File;
use std::path;

use gb_io::reader::SeqReader;
use gb_io::seq::{Feature, QualifierKey};

use failure::Error;

use assembly::*;

fn main() {
    run().unwrap();
}

fn run() -> Result<(), Error> {
    let env = env_logger::Env::default().filter_or(env_logger::DEFAULT_FILTER_ENV, "warn");
    env_logger::Builder::from_env(env).init();
    let args = clap_app!(assembler =>
        (about: "Perform Gibson assembly")
        (version: env!("CARGO_PKG_VERSION"))
        (@arg outdir: +required -o +takes_value "Output directory (created if does not exist)")
        (@arg linear: -l "Include linear products")
        (@arg annotate: -a "Annotate homology regions (don't assemble)")
        (@arg min_homology: -m --("min-homology") default_value("15") "Minimum homology required")
        (@arg FILES: +required ... "Genbank files to use as input")
    )
    .get_matches();
    let files: Vec<_> = values_t_or_exit!(args.values_of("FILES"), String);
    let output = value_t_or_exit!(args.value_of("outdir"), String);
    let min_homology = value_t_or_exit!(args.value_of("min_homology"), usize);
    let names: Vec<_> = files
        .iter()
        .map(|f| path::Path::new(f).file_name())
        .collect();
    let mut seqs = Vec::new();
    let mut seq_names = Vec::new();
    for (fname, name) in files.iter().zip(names) {
        let f = File::open(fname)?;
        for r in SeqReader::new(f) {
            seqs.push(r?);
            seq_names.push(name.unwrap().to_str().unwrap());
        }
    }
    let matches = assembly::find_homology(&seqs, min_homology);
    if args.is_present("annotate") {
        let _ = fs::create_dir_all(&output);
        for (i, mut s) in seqs.into_iter().enumerate() {
            for (&other, positions) in &matches[assembly::Idx(i as u32)] {
                for &Match(start, start_other, len) in positions {
                    let f = Feature {
                        kind: feature_kind!("misc_feature"),
                        location: s.range_to_location(i64::from(start), i64::from(start + len)),
                        qualifiers: vec![(
                            QualifierKey::from("homology_region"),
                            Some(match other {
                                Idx(other) => format!(
                                    "Matches sequence {} position {}",
                                    seq_names[other as usize],
                                    start_other + 1
                                ),
                                IdxRc(other) => format!(
                                    "Matches reverse complement of sequence {} position {}",
                                    seq_names[other as usize],
                                    start_other + 1
                                ),
                            }),
                        )],
                    };
                    s.features.push(f);
                }
            }
            s.write(File::create(
                path::PathBuf::from(&output).join(&seq_names[i]),
            )?)?;
        }
        return Ok(());
    }
    let products = assembly::find_products(&matches, &seqs)?;
    let write = |products: &[assembly::Path], name: Option<&str>| -> Result<(), Error> {
        let path = match name {
            Some(name) => path::PathBuf::from(&output).join(name),
            None => path::PathBuf::from(&output),
        };
        let _ = fs::create_dir_all(&path);
        for (i, p) in products.iter().enumerate() {
            let r = assembly::extract_product_seq(p, &seqs);
            r.write(File::create(path.join(format!("{}.gb", i)))?)?;
        }
        Ok(())
    };
    if args.is_present("linear") {
        write(&products.linear, Some("linear"))?;
        write(&products.circular, Some("circular"))?;
    } else {
        write(&products.circular, None)?;
    }
    Ok(())
}
