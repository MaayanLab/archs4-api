# archs4-api

An API for serving ARCHS4 Data from H5 over S3.

## Technical Details
supervisord launches & monitors aiohttp processes, and nginx process which load balances the aiohttp processes.

## Preparation
For best performance, the h5 file should be created and chunked in a specific way.
- metadata should be created before data tables, in the order that the application wants them:
- meta/genes/gene_symbol
  - defaults might be fine, otherwise large chunks will probably speed things up
- meta/samples/geo_accession
  - defaults might be fine, otherwise large chunks will probably speed things up
- meta/samples/series_id
  - defaults might be fine, otherwise large chunks will probably speed things up
- data/expression
  - chunked using `(n_genes, avg_samples_per_request*2)`
    assuming a shape of `(n_genes, n_samples)`
  - this allows an entire sample to appear in a single chunk. with a chunk size of `avg_samples_per_request*2` each query will take one request to the backend on average (we fetch chunks individually).
- anything else
