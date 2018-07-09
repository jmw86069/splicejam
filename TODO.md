## TODO for jambio

### Portability

* Verify handling of gzip files across architectures, specifically for
makeTx2geneFromGtf(). Fallback plan is to exit gracefully with an informative
message. It appears that data.table::fread() can accept a command as input,
in which case it runs the command and imports the data stream.

### Example data

* Define a small subset of transcript gene models to use in examples
for each function for improved documentation.


### Testing

* include testthat into the package development workflow to help
test and confirm specific input and output criteria.
