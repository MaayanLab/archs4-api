import json
import asyncio
import numpy as np
import pandas as pd
import logging
from pathlib import Path
from aiohttp import web


## Utils ##

bytes_decode = np.vectorize(bytes.decode)

def ensure_list(L):
  if type(L) == list: return L
  else: return [L]

def load_h5(uri):
  import h5py
  if uri.startswith('s3://'):
    import s3fs
    s3 = s3fs.S3FileSystem(anon=True)
    return h5py.File(s3.open(uri, 'rb'))
  else:
    return h5py.File(uri)

async def s3_ctx(app):
  logger = logging.getLogger('s3_ctx')
  logger.info('Setting up h5py over s3...')
  df = load_h5(app['config']['matrix'])
  app['data'] = type('data', tuple(), dict(
    expr=df['data']['expression'],
    genes=bytes_decode(df['meta']['genes']['gene_symbol']),
    series=bytes_decode(df['meta']['samples']['series_id']),
    geo_accession=bytes_decode(df['meta']['samples']['geo_accession']),
  ))
  logger.info('Ready.')
  yield

def serve_mimetype(df, accept=''):
  if accept.startswith('text/tab-separated-values'):
    return web.Response(text=df.to_csv(sep='\t'), content_type='text/tab-separated-values')
  else:
    return web.Response(text=df.to_json(), content_type='application/json')

## ROUTES ##

import yaml
import rororo

operations = rororo.OperationTableDef()

@operations.register('openapi')
@operations.register('openapi.json')
async def openapi_json(request):
  ctx = rororo.get_openapi_context(request)
  return web.json_response(rororo.get_openapi_schema(request.app))

@operations.register('openapi.yaml')
async def openapi_yaml(request):
  ctx = rororo.get_openapi_context(request)
  return web.Response(text=yaml.dump(rororo.get_openapi_schema(request.app)))

@operations.register
async def fetch_data_expression(request):
  ctx = rororo.get_openapi_context(request)
  data = request.app['data']
  #
  if ctx.data.get('genes'):
    gene_filter, = np.where(np.in1d(data.genes, ensure_list(ctx.data['genes'])))
  else:
    gene_filter = slice(None)
  #
  geo_accession = ensure_list(ctx.data['geo_accession'])
  assert len(geo_accession) <= 100, 'Max 100 samples can be processed at a time'
  sample_filter, = np.where(np.in1d(data.geo_accession, geo_accession))
  #
  if isinstance(gene_filter, np.ndarray) and gene_filter.size == 0:
    raise web.HTTPNotFound(reason='No genes provided found in ARCHS4')
  if sample_filter.size == 0:
    raise web.HTTPNotFound(reason='No samples provided found in ARCHS4')
  #
  ret = pd.DataFrame(
    data.expr[:, sample_filter][gene_filter],
    index=data.genes[gene_filter],
    columns=data.geo_accession[sample_filter],
  )
  return serve_mimetype(ret, accept=request.headers.get('Accept', 'application/json'))

@operations.register
async def fetch_meta_genes_gene_symbol(request):
  ctx = rororo.get_openapi_context(request)
  data = request.app['data']
  #
  skip, limit = ctx.parameters.query['skip'], ctx.parameters.query['limit']
  ret = list(data.genes[skip:skip+limit])
  return web.json_response(ret)

@operations.register
async def fetch_meta_samples_geo_accession(request):
  ctx = rororo.get_openapi_context(request)
  data = request.app['data']
  #
  skip, limit = ctx.parameters.query['skip'], ctx.parameters.query['limit']
  ret = list(data.geo_accession[skip:skip+limit])
  return web.json_response(ret)


## CLI ##

import click

@click.command()
@click.option(
  '-m', '--matrix', envvar='ARCHS4_MATRIX',
  default='s3://maayanlab-public/archs4/human_matrix_v10.h5',
  help='ARCHS4 Matrix S3 URI',
)
@click.option(
  '-l', '--listen', envvar='ARCHS4_LISTEN',
  default='5000',
  help='Port or unix socket to listen on'
)
@click.option(
  '-s', '--server-url', envvar='ARCHS4_SERVER_URL',
  default=None,
  help='The server_url'
)
@click.option(
  '-v', '--verbose', envvar='ARCHS4_VERBOSE',
  count=True, default=0,
  help='How verbose this should be, more -v = more verbose'
)
def cli(matrix=None, listen=5000, server_url=None, verbose=0):
  logging.basicConfig(level=30 - (verbose*10))
  app = web.Application()
  app['config'] = dict(matrix=matrix)
  app = rororo.setup_openapi(app,
    Path(__file__).parent/'openapi.yaml',
    operations,
    server_url=server_url,
    use_error_middleware=False,
  )
  app.cleanup_ctx.append(s3_ctx)
  #
  try: listen = json.loads(listen)
  except: pass
  if type(listen) == int:
    web.run_app(app, port=listen)
  elif type(listen) == str:
    web.run_app(app, path=listen)

if __name__ == '__main__':
  from dotenv import load_dotenv
  load_dotenv()
  cli()
