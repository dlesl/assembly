# assembly
This is a library for simulating DNA assembly methods such as [Gibson
assembly](https://en.wikipedia.org/wiki/Gibson_assembly). The
`find_homology` function returns a `MatchArray` representing regions
of sequence homology between any two sequences. The `find_products`
function takes a `MatchArray` and returns any possible linear or
circular products that might be formed by assembling the
sequences. Regions of sequence homology anywhere in the sequence, not
just at its ends, are taken into account. The `extract_product_*`
functions return the predicted sequence.

For more details see the CLI app in `/examples`.

The main aim of this project was to get better at Rust, it hasn't been
at all well tested and probably contains serious bugs! It was largely
inspired (and works similarly to) [pydna](https://github.com/BjornFJohansson/pydna).

## Use it online
* [Simple assembly app](https://dlesl.github.io/clonifier/simple_assembly.html)
* [Clonifier](https://dlesl.github.io/clonifier/)
