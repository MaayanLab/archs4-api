#!/bin/python3

# python client.py series-samples GSE53655 | python client.py fetch > GSE53655.tsv

import re
import json
import click
import urllib.request, urllib.parse, urllib.error
from tqdm import tqdm

base_url = 'https://api.archs4.maayanlab.cloud'
chunk_size = 8192

@click.group()
def cli():
  pass

@cli.command()
@click.option('-o', '--output-file', type=click.File('w'), default='-')
@click.argument('series_id')
def series_samples(series_id, output_file):
  with tqdm(unit=' samples') as t:
    t.set_description('Resolving samples')
    req = urllib.request.Request(f"{base_url}/meta/samples/geo_accession?{urllib.parse.urlencode(dict(series_id=series_id, skip=0, limit=100))}")
    try:
      res = urllib.request.urlopen(req)
    except urllib.error.HTTPError as e:
      raise click.ClickException(message=f"Series ID ({series_id}) Not Found") from e
    start, stop, count = map(int, re.match(r'^(\d+)-(\d+)/(\d+)$', res.info()['content-range']).groups())
    samples = json.load(res)
    print(*samples, sep='\n', file=output_file)
    t.reset(count)
    t.update(stop - start)
    while stop != count:
      req = urllib.request.Request(f"{base_url}/meta/samples/geo_accession?{urllib.parse.urlencode(dict(series_id=series_id, skip=stop, limit=100))}")
      try:
        res = urllib.request.urlopen(req)
      except urllib.error.HTTPError as e:
        raise click.ClickException(message=f"Series ID ({series_id}) Not Found") from e
      start, stop, count = map(int, re.match(r'^(\d+)-(\d+)/(\d+)$', res.info()['content-range']).groups())
      samples = json.load(res)
      print(*samples, sep='\n', file=output_file)
      t.update(stop - start)

@cli.command()
@click.option('-i', '--samples-from', type=click.File('r'), default='-')
@click.option('-o', '--output-file', type=click.File('w'), default='-')
def fetch(samples_from, output_file):
  with tqdm(unit=' sample') as t:
    t.set_description('Loading samples')
    geo_accession = [
      sample_id
      for sample_id in map(str.strip, samples_from)
      if sample_id and t.update(1) is None
    ]
    t.reset(total=len(geo_accession) + 1)
    t.set_description('Fetching samples')
    req = urllib.request.Request(
      f"{base_url}/data/expression/T",
      headers={
        'Accept': 'text/tsv',
        'Content-Type': 'application/json',
      },
      data=json.dumps(dict(
        geo_accession=geo_accession,
      )).encode('utf-8'),
    )
    res = urllib.request.urlopen(req)
    while chunk := res.read(chunk_size):
      chunk = chunk.decode('utf-8')
      output_file.write(chunk)
      t.update(chunk.count('\n'))

if __name__ == '__main__':
  cli()
