# kbo-cli
Command-line interface for the [kbo](https://github.com/tmaklin/kbo) local aligner.

Documentation is available at [https://docs.rs/kbo](https://docs.rs/kbo).

## Installation
### Download
Binaries are available from the [Releases page](https://github.com/tmaklin/kbo-cli/releases).

### From bioconda
Run
``` text
conda install -c bioconda -y kbo-cli
```

### Using cargo
Run `cargo install kbo-cli`.

### Compile from source
Run
``` text
git clone https://github.com/tmaklin/kbo-cli
cd kbo-cli
cargo build --release
```
This will build the `kbo` executable in `target/release/kbo` directory.

## Usage
kbo-cli provides access to three main operations:


- `kbo call` calls single and multi base substitutions,
  insertions, and deletions in a query sequence against a reference and
  reports their positions and sequences. Call is useful for problems that
  require [.vcf files](https://samtools.github.io/hts-specs/VCFv4.2.pdf).
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
