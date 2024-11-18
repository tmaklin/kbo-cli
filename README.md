# kbo-cli
Command-line interface for the [kbo](https://docs.rs/kbo) local aligner.

Documentation is available at [https://docs.rs/kbo](https://docs.rs/kbo).

## Installation
### Cargo
Run `cargo install kbo-cli`.

### From source
Run
``` text
git clone https://github.com/tmaklin/kbo-cli
cd kbo-cli
cargo build --release
```
This will build the `kbo` executable in `target/release/kbo` directory.

## Usage
kbo-cli provides access to two main operations:

- `kbo find` matches the _k_-mers in a query sequence with the
  reference and reports the local alignment segments found within the
  reference. Find is useful for problems that can be solved with
  [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi).
- `kbo map` maps the query sequence against a reference
  sequence, and reports the nucleotide sequence of the alignment relative to
  the reference. Map solves the same problem as
  [snippy](https://github.com/tseemann/snippy) and [ska
  map](https://docs.rs/ska/latest/ska/#ska-map).

For usage instructions, see the documentation at [https://docs.rs/kbo](https://docs.rs/kbo).

## License
kbo-cli is dual-licensed under the [MIT](LICENSE-MIT) and [Apache 2.0](LICENSE-APACHE) licenses.
