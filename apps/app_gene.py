import pickle, copy, os, datetime, subprocess, json
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy.stats import entropy
import time
from io import StringIO

import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_table_experiments as dt
from dash.dependencies import Input, Output, State
import plotly.graph_objs as go
import flask
import plotly
from flask_caching import Cache

import inDelphi
import generalStats
import lib, header

import boto3, botocore
import os
s3 = boto3.resource('s3', aws_access_key_id = os.environ['S3_KEY'], aws_secret_access_key = os.environ['S3_SECRET'])

from indelphi_app import app

# init
if not os.path.isdir('local-s3/'):
  os.mkdir('local-s3/')
else:
  subprocess.check_output('rm -rf local-s3/*', shell = True)

# Set up flask caching
CACHE_CONFIG = {
  'CACHE_TYPE': 'redis',
  'CACHE_REDIS_URL': os.environ.get('REDIS_URL', '')
}
cache = Cache()
cache.init_app(app.server, config = CACHE_CONFIG)
cache_timeout = 120

# Remove these plotly modebar buttons to limit interactivity
modebarbuttons_2d = ['zoom2d', 'pan2d', 'select2d', 'lasso2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'resetScale2d', 'hoverClosestCartesian', 'hoverCompareCartesian', 'toggleSpikelines']

## Parameters

###################################################################
###################################################################
##
# App layout
##
layout = html.Div([

  ###################################################
  # Hidden divs for light data storage
  ###################################################
  html.Div(
    [
      html.Div(
        id = 'G_hidden-pred-df-stats-signal',
        children = 'init'
      ),
      html.Div(
        id = 'G_table-stats-signal',
        children = 'init'
      ),
      html.Div(
        id = 'G_hidden-selected-genome',
        children = 'init'
      ),
      html.Div(
        id = 'G_hidden-selected-gene',
        children = 'init'
      ),
      html.Div(
        id = 'G_hidden-cache-submit-button',
        children = '%s' % (time.time())
      ),
      html.Div(
        id = 'G_hidden-sort-module-interaction',
        children = '%s' % (time.time())
      ),
      html.Div(
        id = 'G_hidden-clickData',
        children = '%s init' % (time.time())
      ),
      html.Div(
        id = 'G_hidden-selected-id',
        children = ''
      ),

      # Datatable
      dt.DataTable(
        id = 'G_table-stats',
        rows = [{}], # init rows
        selected_row_indices = [],
      ),

      dcc.Location(
        id = 'G_url',
        refresh = False,
      ),
    ],
    style = dict(
      display = 'none',
    ),
  ),

  ###################################################
  # Header
  ###################################################
  html.Div(
    [
      ###################################################
      # Upper header
      ###################################################
      header.get_navigation_header('gene'),

      ###################################################
      # Genome choice
      ###################################################
      html.Div(
        [
          html.Div(
            [
              # Left
              html.Div(
                [
                  html.Span('Genome: '),
                ],
                style = dict(
                  display = 'table-cell',
                  textAlign = 'right',
                  width = '50%',
                  transform = 'translateX(-10px)',
                ),
              ),
              # Middle
              html.Div(
                [
                  dcc.RadioItems(
                    id = 'G_genome-radio',
                    options = [
                      {'label': 'Human (hg38)', 'value': 'hg38'},
                      {'label': 'Mouse (mm10)', 'value': 'mm10'},
                    ],
                    value = 'hg38'
                  )
                ],
                style = dict(
                  display = 'table-cell',
                  width = '30%',
                ),
              ),
              # Right
              html.Div(
                [],
                style = dict(
                  display = 'table-cell',
                  textAlign = 'left',
                  width = '20%',
                  transform = 'translateX(10px)',
                ),
              ),
            ],
            style = dict(
              display = 'table-row',
            ),
          ),
        ],
        style = dict(
          display = 'table',
          width = '100%',
          marginBottom = 10,
        ),
      ),

      ###################################################
      # Gene dropdown
      ###################################################
      html.Div(
        [
          html.Div(
            [
              # Left
              html.Div(
                [
                  html.Span('Gene: '),
                ],
                style = dict(
                  display = 'table-cell',
                  textAlign = 'right',
                  width = '50%',
                  transform = 'translateX(-10px)',
                ),
              ),
              # Middle
              html.Div(
                [
                  dcc.Dropdown(
                    id = 'G_gene-dropdown',
                    placeholder = 'Type to search for a gene',
                  ),
                ],
                style = dict(
                  display = 'table-cell',
                  width = '25%',
                ),
              ),
              # Right
              html.Div(
                [],
                style = dict(
                  display = 'table-cell',
                  textAlign = 'left',
                  width = '25%',
                  transform = 'translateX(10px)',
                ),
              ),
            ],
            style = dict(
              display = 'table-row',
            ),
          ),
        ],
        style = dict(
          display = 'table',
          width = '100%',
          marginBottom = 10,
        ),
      ),

      ###################################################
      # Cell type
      ###################################################
      html.Div(
        [
          html.Div(
            [
              # Left
              html.Div(
                [
                  html.Span('Cell type: '),
                ],
                style = dict(
                  display = 'table-cell',
                  textAlign = 'right',
                  width = '50%',
                  transform = 'translateX(-10px)',
                ),
              ),
              # Middle
              html.Div(
                [
                  dcc.Dropdown(
                    options = [
                      {'label': 'HCT116', 'value': 'HCT116', 'disabled': True},
                      {'label': 'HEK293', 'value': 'HEK293', 'disabled': True},
                      {'label': 'K562', 'value': 'K562', 'disabled': True},
                      {'label': 'mESC', 'value': 'mESC'},
                      {'label': 'U2OS', 'value': 'U2OS'},
                    ],
                    id = 'G_celltype_dropdown',
                    searchable = False,
                    clearable = False,
                    value = 'mESC',
                  ),
                ],
                style = dict(
                  display = 'table-cell',
                  width = '10%',
                ),
              ),
              # Right
              html.Div(
                [
                  html.Div(
                    [
                      html.Img(
                        src = '/staticfiles/tooltip_logo',
                        className = 'tooltiprightlogo',
                      ),
                      html.Span(
                        'Choose a cell type specific version of inDelphi. If your cell type of interest is not listed here, we recommend using mESC if your cell type has no expected defects in DNA repair. Contradicting the genome choice is not a problem: for example, human embryonic stem cells are likely to have more similar DNA repair outcomes to mESC than human cancer cell lines.',
                        className = 'tooltiprighttext',
                        style = dict(width = '200px',)
                      ),
                    ], 
                    className = 'tooltipright',
                  ),
                ],
                style = dict(
                  display = 'table-cell',
                  textAlign = 'left',
                  width = '40%',
                  transform = 'translateX(10px)',
                ),
              ),
            ],
            style = dict(
              display = 'table-row',
            ),
          ),
        ],
        style = dict(
          display = 'table',
          width = '100%',
          marginBottom = 10,
        ),
      ),


      ###################################################
      # Submit button
      ###################################################
      # Submit button
      html.Div([
        html.Button(
          'SUBMIT',
          id = 'G_submit_button',
          style = dict(
            boxShadow = '1px 3px 6px 0 rgba(0, 0, 0, 0.2)',
          ),
        )],
        style = dict(
          textAlign = 'center',
          marginBottom = 15,
        ),
      ),

    ],
    style = dict(
      backgroundColor = 'white',
      width = '1010px',
      position = 'relative',
      left = '50%',
      transform = 'translate(-50%, 0px)',
      borderBottom = '3px solid #777777',
      marginBottom = '50px',
    ),
  ),

  ###################################################
  # Post-computation settings module + Histograms (sticky)
  ###################################################
  html.Div(
    [
      # Module
      html.Div([

        # Header
        html.Div([
          html.Div([
            html.Strong('',
              id = 'G_postcomp_module_header',
            )],
            className = 'module_header_text'),
          ],
          className = 'module_header'
        ),

        # Module body
        html.Div(
          [
            # Row: Display kgIDs...
            html.Div(
              [
                html.Strong(
                  'Display kgIDs:',
                  style = dict(
                    textAlign = 'right',
                    marginRight = '5px',
                    height = '36px',  # height of one dropdown line
                    lineHeight = '36px',  # centers vertically
                  ),
                  className = 'three columns',
                ),

                # Multi drop down to select columns
                dcc.Dropdown(
                  id = 'G_dropdown-kgid',
                  multi = True,
                  searchable = False,
                  clearable = False,
                  className = 'nine columns',
                ),
              ],
              style = dict(
                # width = '1050px',
                marginBottom = '5px',
                marginTop = '10px',
              ),
              className = 'row',
              id = 'G_row_dropdown-kgid',
            ),

            # Row: Display columns...
            html.Div(
              [
                html.Strong(
                  'Display columns:',
                  style = dict(
                    textAlign = 'right',
                    marginRight = '5px',
                    height = '36px',  # height of one dropdown line
                    lineHeight = '36px',  # centers vertically
                  ),
                  className = 'three columns',
                ),

                # Multi drop down to select columns
                dcc.Dropdown(
                  id = 'G_dropdown-columns',
                  options = [
                    {'label': 'Exon number', 'value': 'Exon number'},
                    {'label': 'Distance to 5\' exon boundary', 'value': 'Dist. to 5\' end'},
                    {'label': 'Distance to 3\' exon boundary', 'value': 'Dist. to 3\' end'},
                    {'label': 'Precision', 'value': 'Precision'},
                    {'label': 'Frameshift (%)', 'value': 'Frameshift (%)'},
                    {'label': 'Frame +0 (%)', 'value': 'Frame +0 (%)'},
                    {'label': 'Frame +1 (%)', 'value': 'Frame +1 (%)'},
                    {'label': 'Frame +2 (%)', 'value': 'Frame +2 (%)'},
                    {'label': 'Microhomology strength', 'value': 'MH strength'},
                    {'label': 'Most frequent genotype (%)', 'value': 'M.F. gt (%)'},
                    {'label': 'Most frequent deletion (%)', 'value': 'M.F. del (%)'},
                    {'label': 'Most frequent insertion (%)', 'value': 'M.F. ins (%)'},
                    {'label': 'Expected indel length', 'value': 'Exp. indel len'},
                  ],
                  multi = True,
                  searchable = False,
                  clearable = False,
                  value = ['Exon number', 'Dist. to 5\' end', 'Dist. to 3\' end', 'Precision', 'Frameshift (%)', 'Frame +0 (%)'],
                  className = 'nine columns',
                ),
              ],
              style = dict(
                # width = '1050px',
                marginBottom = '5px',
                marginTop = '10px',
              ),
              className = 'row',
              id = 'G_row_dropdown-columns',
            ),

            # Row: Sort by...
            html.Div(
              [
                html.Strong(
                  'Sort by:  ',
                  className = 'three columns',
                  style = dict(
                    textAlign = 'right',
                    marginRight = '5px',
                    height = '36px',
                    lineHeight = '36px',
                  ),
                ),
                # Sorting columns
                dcc.Dropdown(
                  id = 'G_dropdown-sortcol',
                  options = [],
                  searchable = False,
                  clearable = False,
                  className = 'three columns',
                ),
                # Sort direction
                dcc.RadioItems(
                  id = 'G_sortdirection',
                  options = [
                    {'label': 'Ascending', 'value': 'Ascending'},
                    {'label': 'Descending', 'value': 'Descending'},
                  ],
                  value = 'Descending',
                  labelStyle = {'display': 'inline-block'},
                  className = 'six columns',
                  style = dict(
                    marginLeft = 5,
                    height = '36px',
                    lineHeight = '36px',
                  ),
                ),
              ],
              style = dict(
                marginBottom = '10px',
              ),
              className = 'row',
              id = 'G_row_dropdown-sortcol',
            ),

            # Links
            html.Div([
              html.Div(
                # Sharable link
                html.A(
                  '🔗 Shareable link to page before computation', 
                  id = 'G_page-link'
                )
              ),
              html.Div(
                # Download link: summary statistics
                html.A(
                  '📑 Download table of predictions',
                  id = 'G_download-link'
                )
              ),
              html.Div([
                html.Span(
                  'Note: Online visualization is limited to 1000 gRNAs.',
                )
              ])
            ], style = dict(
                textAlign = 'center',
                height = 90,
              )
            ),

          ],
        ),

        ##
        ], 
        style = dict(
          transform = 'translateX(240px)',
          width = '970px',
          boxShadow = '1px 3px 6px 0 rgba(0, 0, 0, 0.2)',
          marginBottom = '50px',
          position = 'relative',
          zIndex = 10,
        ),
      ),

      # Hists
      html.Div(
        dcc.Graph(
          id = 'G_hist-stats',
          config = dict(
            modeBarButtonsToRemove = modebarbuttons_2d,
            displaylogo = False,
            displayModeBar = False,
          ),
        ),
        id = 'G_hist-stats-div',
        style = dict(
          display = 'none',
          position = 'relative',
          zIndex = 1,
        )
      ),
    ],
    # body style
    id = 'G_postcomputation_settings',
    className = 'batch_postcomputation_sticky',
    style = dict(
      display = 'none',
    ),
  ),

  ###################################################
  # Plots
  ###################################################
  html.Div(
    [
      # Plots
      html.Div(
        dcc.Graph(
          id = 'G_plot-stats',
          config = dict(
            modeBarButtonsToRemove = modebarbuttons_2d,
            displaylogo = False,
            displayModeBar = False,
          ),
        ),
        id = 'G_plot-stats-div',
        style = dict(
          display = 'none',
        ),
        className = 'animate-bottom',
      ),

    ],
    # body style
    style = dict(
    ),
  ),
  ##

  ],  # body div
  style = dict(
    # width = '1150px',
    width = '1450px',
    margin = '0 auto',
  )
)

#######################################################################
#########################      CALLBACKS      #########################
#######################################################################

##
# Hidden button callbacks
##
@app.callback(
  Output('G_hidden-cache-submit-button', 'children'),
  [Input('G_submit_button', 'n_clicks')])
def update_submit_button_time(n_clicks):
  return '%s' % (time.time())

@app.callback(
  Output('G_hidden-sort-module-interaction', 'children'),
  [Input('G_row_dropdown-columns', 'n_clicks'),
   Input('G_row_dropdown-sortcol', 'n_clicks')])
def update_sort_time(v1, v2):
  return '%s' % (time.time())

@app.callback(
  Output('G_hidden-clickData', 'children'),
  [Input('G_plot-stats', 'clickData')])
def update_hidden_clickdata(clickData):
  return '%s %s' % (time.time(), clickData['points'][0]['pointNumber'])

##
# URL callbacks
##
@app.callback(
  Output('G_genome-radio', 'value'),
  [Input('G_url', 'pathname')],
  [State('G_genome-radio', 'value')])
def update_genome_build_from_url(url, default_value):
  valid_flag, dd = lib.parse_valid_url_path_gene(url)
  if valid_flag:
    return dd['genome_build']
  return default_value

@app.callback(
  Output('G_gene-dropdown', 'value'),
  [Input('G_url', 'pathname')],
  [State('G_gene-dropdown', 'value')])
def update_gene_from_url(url, default_value):
  valid_flag, dd = lib.parse_valid_url_path_gene(url)
  if valid_flag:
    return dd['gene']
  return default_value

@app.callback(
  Output('G_celltype_dropdown', 'value'),
  [Input('G_url', 'pathname')],
  [State('G_celltype_dropdown', 'value')])
def update_celltype_from_url(url, default_value):
  valid_flag, dd = lib.parse_valid_url_path_gene(url)
  if valid_flag:
    return dd['celltype']
  return default_value

@app.callback(
  Output('G_dropdown-sortcol', 'value'),
  [Input('G_dropdown-sortcol', 'options')],
  [State('G_dropdown-sortcol', 'value'),
   State('G_url', 'pathname')])
def update_sortcols_from_url(options, default_value, url):
  valid_flag, dd = lib.parse_valid_url_path_gene(url)
  if not valid_flag or dd['sort_by'] == '-':
    return default_value
  else:
    all_options = [s['value'] for s in options]
    idx = int(dd['sort_by'])
    return sorted(all_options)[idx]

@app.callback(
  Output('G_sortdirection', 'value'),
  [Input('G_url', 'pathname')],
  [State('G_sortdirection', 'value')])
def update_sortdir_from_url(url, default_value):
  valid_flag, dd = lib.parse_valid_url_path_gene(url)
  if valid_flag:
    return dd['sort_dir']
  else:
    return default_value

@app.callback(
  Output('G_dropdown-columns', 'value'),
  [Input('G_url', 'pathname')],
  [State('G_dropdown-columns', 'value'),
   State('G_dropdown-columns', 'options')])
def update_columns_from_url(url, default_value, options):
  all_options = [s['value'] for s in options]
  valid_flag, dd = lib.parse_valid_url_path_gene(url)
  if valid_flag:
    value = []
    alphabetical_options = sorted(all_options)
    for idx, flag in enumerate(dd['chosen_columns']):
      if flag == '1':
        value.append(alphabetical_options[idx])
    return value
  else:
    return default_value


##
# Header callbacks
##
@app.callback(
  Output('G_gene-dropdown', 'options'),
  [Input('G_genome-radio', 'value')])
def update_gene_dropdown_choices(genome_build):
  stats_dir = os.path.dirname(os.path.realpath(__file__)) + '/statistics/'
  if genome_build == 'mm10':
    return generalStats.mm10_choices
  elif genome_build == 'hg38':
    return generalStats.hg38_choices

@app.callback(
  Output('G_submit_button', 'children'),
  [Input('G_gene-dropdown', 'value')],
  [State('G_submit_button', 'children')])
def update_submit_button_text(selected_gene, prev_value):
  if selected_gene is None:
    return 'SELECT A GENE'
  else:
    return 'SUBMIT'


@app.callback(
  Output('G_submit_button', 'style'),
  [Input('G_gene-dropdown', 'value')],
  [State('G_submit_button', 'style')])
def update_submit_button_style(selected_gene, style):
  if selected_gene is None:
    style['backgroundColor'] = '#86898C'
    style['color'] = 'white'
  else:
    style['backgroundColor'] = '#00A0DC'
    style['color'] = 'white'
  return style


##
# AWS S3 download callback
##
@cache.memoize()
def grab_s3_stats_cache(parameters):
  genome_build, gene, celltype = parameters

  query_fn = '%s_%s_SpCas9_%s.csv' % (genome_build, celltype, gene)
  local_dir = 'local-s3/'
  s3.Bucket('indelphi-storage').download_file(query_fn, local_dir + query_fn)

  all_stats = pd.read_csv(local_dir + query_fn, index_col = 0)
  all_stats['ID'] = all_stats.index + 1
  all_stats['PAM'] = [s[63:66] for s in all_stats['Local context']]
  all_stats['MH strength'] = np.log(all_stats['Phi'])

  dd = defaultdict(list)
  for idx, row in all_stats.iterrows():
    sm_link = lib.encode_dna_to_url_path_single(row['Local context'], 60, celltype)
    dd['URL'].append('%s' % (sm_link))

    if row['Exon strand'] == row['gRNA strand w.r.t. exon strand']:
      dd['Strand'].append('+')
    else:
      dd['Strand'].append('-')

    if row['Exon strand'] == '+':
      cutsite_coord = int(row['Exon start']) + int(row['Cutsite distance to 5p boundary'])
    else:
      # for col in all_stats.columns:
        # print(col, row[col])
      cutsite_coord = int(row['Exon start']) + int(row['Cutsite distance to 3p boundary'])
    dd['Cutsite coordinate'].append(cutsite_coord)

  for col in dd:
    all_stats[col] = dd[col]

  all_stats['Distance to 5\' exon boundary'] = all_stats['Cutsite distance to 5p boundary']
  all_stats['Distance to 3\' exon boundary'] = all_stats['Cutsite distance to 3p boundary']

  return all_stats

@app.callback(
  Output('G_hidden-pred-df-stats-signal', 'children'),
  [Input('G_submit_button', 'n_clicks')],
  [State('G_genome-radio', 'value'),
   State('G_gene-dropdown', 'value'),
   State('G_celltype_dropdown', 'value')]
)
def update_df_stats(n_clicks, genome_build, gene, celltype):
  parameters = (genome_build, gene, celltype)
  grab_s3_stats_cache(parameters)
  return parameters

##
# Module header callbacks, Advanced options hiding/showing
##
@app.callback(
  Output('G_hidden-selected-genome', 'children'),
  [Input('G_table-stats-signal', 'children')],
  [State('G_genome-radio', 'value')]
)
def update_hidden_selected_genome(signal, genome):
  return genome

@app.callback(
  Output('G_hidden-selected-gene', 'children'),
  [Input('G_table-stats-signal', 'children')],
  [State('G_gene-dropdown', 'value')]
)
def update_hidden_selected_gene(signal, gene):
  return gene

@app.callback(
  Output('G_postcomp_module_header', 'children'),
  [Input('G_table-stats-signal', 'children'),
   Input('G_hidden-selected-genome', 'children'),
   Input('G_hidden-selected-gene', 'children')]
)
def update_postcomp_module_header(table_signal, genome_build, gene):
  df = make_table_stats_cache(table_signal)
  return 'Results of %s SpCas9 (NGG) gRNAs targeting %s in %s' % (len(df), gene, genome_build)

##
# kgID, column selection and sorting callbacks
##
@app.callback(
  Output('G_dropdown-sortcol', 'options'),
  [Input('G_dropdown-columns', 'value')])
def update_sortcol_options(values):
  options = []
  for value in values:
    options.append({'label': value, 'value': value})
  return options

@app.callback(
  Output('G_dropdown-kgid', 'options'),
  [Input('G_dropdown-kgid', 'value')],
  [State('G_hidden-pred-df-stats-signal', 'children')]
)
def update_dropdown_kgid_options(value, signal):
  if signal == 'init':
    assert False, 'init'
  stats = grab_s3_stats_cache(signal)
  kgids = list(set(stats['kgID']))
  sizes = [len(stats[stats['kgID'] == kgid]) for kgid in kgids]
  options = []
  total_size_of_selected = sum([sizes[kgids.index(s)] for s in value])
  for kgid, size in zip(kgids, sizes):
    curr_opt = {'label': '%s (%s gRNAs)' % (kgid, size), 'value': kgid}
    if kgid not in value:
      if size + total_size_of_selected > 1000:
        curr_opt['disabled'] = True
    options.append(curr_opt)
  return options

@app.callback(
  Output('G_dropdown-kgid', 'value'),
  [Input('G_hidden-pred-df-stats-signal', 'children')]
)
def update_dropdown_kgid_value(signal):
  if signal == 'init':
    assert False, 'init'
  stats = grab_s3_stats_cache(signal)
  kgids = set(stats['kgID'])
  sizes = [len(stats[stats['kgID'] == kgid]) for kgid in kgids]
  kgids_sorted = [x for _,x in sorted(zip(sizes, kgids), reverse = True)]
  sizes_sorted = sorted(sizes, reverse = True)

  # Select the largest possible
  for idx in range(len(sizes_sorted)):
    if sizes_sorted[idx] > 1000:
      sizes_sorted = sizes_sorted[1:]
      kgids_sorted = kgids_sorted[1:]
    else:
      break

  for idx in range(1, len(sizes_sorted)):
    if sum(sizes_sorted[:idx]) > 1000:
      return kgids_sorted[:idx - 1]
  return kgids_sorted

##
# Stats table callbacks
## 
@cache.memoize(timeout = cache_timeout)
def make_table_stats_cache(parameters):
  parameters = json.loads(parameters)
  signal, chosen_columns, sort_col, sort_direction, kgids = parameters

  stats = grab_s3_stats_cache(signal)

  # Drop unselected kgids
  stats = stats[stats['kgID'].isin(kgids)]
  assert len(stats) <= 1000

  # Drop extra cols
  drop_cols = [
    '1-bp ins frequency',
    'MH del frequency',
    'MHless del frequency',
  ]
  stats = stats.drop(drop_cols, axis = 1)

  # Rename to shorter versions
  stats = lib.rename_batch_columns(stats)

  # Sort by, if possible
  if sort_col is not None and sort_direction is not None:
    if sort_direction == 'Ascending':
      ascending_flag = True
    else:
      ascending_flag = False
    stats = stats.sort_values(by = sort_col, ascending = ascending_flag)

  # Reformat floats
  stats_cols = list(stats.columns)
  nonstat_cols = [
    'ID',
    'PAM',
    'URL',
    'Cas9 type',
    'Celltype',
    'Chromosome',
    'Cutsite distance to 3p boundary',
    'Cutsite distance to 5p boundary',
    'Exon end',
    # 'Exon number',
    'Exon start',
    'Exon strand',
    'Gene symbol',
    'Genome',
    'Local context',
    'Local cutsite',
    'gRNA',
    'gRNA strand w.r.t. exon strand',
    'kgID',
    'Strand',
    'Cutsite coordinate',
  ]

  for nonstat_col in nonstat_cols:
    stats_cols.remove(nonstat_col)
  for stat_col in stats_cols:
    # Filter down to selected columns
    if stat_col not in chosen_columns:
      stats.drop(stat_col, axis = 1, inplace = True)
      continue
    # Reformat
    if stat_col in ['Precision', 'MH strength']:
      stats[stat_col] = [float('%.2f' % (s)) for s in stats[stat_col]]    
    else:
      stats[stat_col] = [float('%.1f' % (s)) for s in stats[stat_col]]    

  # Reorder columns
  stats = stats[nonstat_cols + lib.order_chosen_columns(chosen_columns)]
  stats = stats.reset_index(drop = True)
  return stats

@app.callback(
  Output('G_table-stats-signal', 'children'), 
  [Input('G_hidden-pred-df-stats-signal', 'children'),
   Input('G_dropdown-columns', 'value'),
   Input('G_dropdown-sortcol', 'value'),
   Input('G_sortdirection', 'value'),
   Input('G_dropdown-kgid', 'value'),
  ])
def update_stats_table(signal, chosen_columns, sort_col, sort_direction, kgids):
  if signal == 'init':
    assert False, 'init'

  parameters = (signal, chosen_columns, sort_col, sort_direction, kgids)
  parameters = json.dumps(parameters)
  make_table_stats_cache(parameters)
  return parameters


@app.callback(
  Output('G_table-stats', 'selected_row_indices'),
  [Input('G_hidden-clickData', 'children'),
   Input('G_hidden-cache-submit-button', 'children'),
   Input('G_dropdown-columns', 'value'),
   Input('G_dropdown-sortcol', 'value'),
   Input('G_table-stats-signal', 'children')],
  [State('G_table-stats', 'selected_row_indices'),
   State('G_hidden-sort-module-interaction', 'children'),
   State('G_hidden-selected-id', 'children'),
   State('G_url', 'pathname'),
   State('G_postcomputation_settings', 'n_clicks'),
   State('G_plot-stats-div', 'n_clicks'),
   State('G_submit_button', 'n_clicks'),
   ])
def update_statstable_selected(clickData, submit_time, col_values, sortcol_value, table_signal, selected_row_indices, sort_time, prev_id, url, nc1, nc2, nc_submit):
  if not bool(nc1 and nc2) and nc_submit == 1:
    # On page load, select row from URL
    valid_flag, dd = lib.parse_valid_url_path_gene(url)
    if valid_flag:
      if dd['row_select'] != '-':
        return [int(dd['row_select'])]

  # Only allow selecting one point in plot-stats
  submit_time = float(submit_time)
  sort_time = float(sort_time)
  click_time = float(clickData.split()[0])
  click_idx = clickData.split()[1]
  if click_idx == 'init':
    return []
  else:
    click_idx = int(click_idx)

  submit_intxn = bool(submit_time > max(sort_time, click_time))
  click_intxn = bool(click_time > max(sort_time, submit_time))
  sort_intxn = bool(sort_time > max(click_time, submit_time))

  print('Submit: %s' % (submit_intxn))
  print('Click: %s' % (click_intxn))
  print('Sort: %s' % (sort_intxn))

  if sort_intxn and prev_id != '':
    # If changing sort col or direction, clear the selected rows. Otherwise, the wrong row is selected after sorting. Preferably, keep the selected row and update the index.
    selected_row_indices = []
    df = make_table_stats_cache(table_signal)

    # new_idx = int(df[df['ID'] == int(prev_id)].index[0])
    id_list = list(df['ID'])
    real_new_idx = id_list.index(int(prev_id))
    display_new_idx = len(df) - real_new_idx - 1
    new_idx = display_new_idx
    selected_row_indices = [new_idx]
  elif submit_intxn:
    # if hitting submit button, clear the selected rows. Otherwise, selecting a row M > number of rows N in new query, will fail
    selected_row_indices = []
  elif click_intxn:
    # Must be triggered by clickData
    # Update selections in table based on clicking plot
    if selected_row_indices != [click_idx]:
      selected_row_indices = [click_idx]
    else:
      # Point already selected, user clicked on same point twice:
      # so, deselect
      selected_row_indices = []
  return selected_row_indices

@app.callback(
  Output('G_hidden-selected-id', 'children'),
  [Input('G_table-stats', 'selected_row_indices')],
  [State('G_table-stats-signal', 'children')])
def update_hidden_selected_id(selected_idx, table_signal):
  if len(selected_idx) == 0:
    return ''
  idx = selected_idx[0]
  df = make_table_stats_cache(table_signal)
  return list(df['ID'])[idx]


##
# Plot stats callback: styles, hide when no figure
##
@app.callback(
  Output('G_plot-stats-div', 'style'),
  [Input('G_plot-stats', 'figure')])
def update_stats_plot_style(fig):
  if fig is None:
    return {'display': 'none'}
  else:
    return {}

@app.callback(
  Output('G_hist-stats-div', 'style'),
  [Input('G_hist-stats', 'figure')])
def update_hist_plot_style(fig):
  if fig is None:
    return {'display': 'none'}
  else:
    return {}

@app.callback(
  Output('G_postcomputation_settings', 'style'),
  [Input('G_plot-stats', 'figure')])
def update_postcomputation_settings_style(fig):
  if fig is None:
    return {'display': 'none'}
  else:
    return {}


########################################################
# Plot stats callback
########################################################
@app.callback(
    Output('G_plot-stats', 'figure'),
    [Input('G_table-stats-signal', 'children'),
     Input('G_table-stats', 'selected_row_indices')])
def update_stats_plot(table_signal, selected_row_indices):
  df = make_table_stats_cache(table_signal)
  # Determine statistics to plot
  stats_cols = lib.order_chosen_columns(list(df.columns))

  fig = plotly.tools.make_subplots(
    rows = 1, cols = len(stats_cols),
    shared_yaxes = True)

  # Color selected markers
  if len(selected_row_indices) > 0:
    selected_row_index = selected_row_indices[0]
  else:
    selected_row_index = None
  selected_line = dict()

  yrange = np.arange(1, len(df.index) + 1)

  # Generate each plot
  for idx, stats_col in enumerate(stats_cols):
    subplot_num = idx + 1
    marker = {'color': [lib.get_color(stats_col)] * len(df)}
    for i in (selected_row_indices or []):
      marker['color'][i] = '#000000'
    # Gray lines
    fig.append_trace(
      go.Bar(
        x = df[stats_col][::-1],
        y = yrange,
        orientation = 'h',
        hoverinfo = 'skip',
        width = 0.1,
        opacity = 0.2,
        marker = dict(
          color = 'gray',
        )
      ), 
      1, subplot_num
    )

    # Scatter
    fig.append_trace(
      go.Scattergl(
        x = df[stats_col][::-1],
        y = yrange,
        mode = 'markers',
        marker = marker,
        name = '',
      ), 
      1, subplot_num
    )

    if selected_row_index is not None:
      selected_line[subplot_num] = (df.index[selected_row_index], df[stats_col][len(df) - selected_row_index - 1])

  # Format y tick texts: ID, gRNA, PAM, orientation, URL.
  yticktexts = []
  fw_ids = lib.get_fixedwidth_ID(df['ID'])
  fw_kgids = lib.get_fixedwidth_items(df['kgID'])
  fw_coords = lib.get_fixedwidth_items(df['Cutsite coordinate'])
  for idx, row in df.iterrows():
    row_text = '%s %s %s %s %s %s <a href="%s">details</a> %s' % (row['gRNA'], row['PAM'], row['Chromosome'], fw_coords[idx], row['Strand'], fw_kgids[idx], row['URL'], fw_ids[idx])
    yticktexts.append(row_text)


  # Subplot formatting
  fig['layout']['barmode'] = 'stack'
  fig['layout']['yaxis1'].update(
    fixedrange = True,
    # autorange = False,
    tickvals = yrange,
    range = [min(yrange) - 1, max(yrange) + 1],
    ticktext = yticktexts[::-1],
    tickfont = dict(
      size = 12,
      family = 'monospace',
    ),
    zeroline = True,
    zerolinewidth = 2,
    # autorange = 'reversed',
  )

  all_shapes = []
  x_domains = lib.get_x_domains(len(stats_cols))
  for idx, stats_col in enumerate(stats_cols):
    subplot_num = idx + 1
    [xmin, xmax] = lib.get_batch_statcol_xrange(df[stats_col], stats_col)
    fig['layout']['xaxis%s' % (subplot_num)].update(
      # title = stats_col,
      domain = x_domains[idx],
      fixedrange = True,
      # showgrid = False,
      showgrid = True,
      zeroline = False,
      titlefont = dict(
        size = 12,
      ),
      range = [xmin, xmax],
      # showspikes = True,
      # spikesnap = 'cursor',
      # spikemode = 'across+marker',
      # spikedash = 'solid',
      # spikethickness = 1,
      # spikecolor = '#777',
    )

    if selected_row_index is not None:
      all_shapes.append(
        lib.get_batch_select_line(
          x0 = selected_line[subplot_num][1],
          x1 = selected_line[subplot_num][1],
          y0 = 0,
          y1 = len(df),
          xref = 'x%s' % (subplot_num),
          yref = 'y1',
        )
      )
      all_shapes.append(
        lib.get_batch_select_line(
          x0 = xmin,
          x1 = xmax,
          y0 = selected_line[subplot_num][0] + 1,
          y1 = selected_line[subplot_num][0] + 1,
          xref = 'x%s' % (subplot_num),
          yref = 'y1',
        )
      )

  fig['layout']['shapes'] = all_shapes
  # Global figure formatting
  fig['layout']['showlegend'] = False
  fig['layout']['hovermode'] = 'y'
  # fig['layout']['spikedistance'] = -1
  fig['layout']['width'] = 455 + len(stats_cols) * 150
  fig['layout']['height'] = 150 + len(df) * 11
  fig['layout']['margin'] = {
    'l': 430,
    'r': 25,
    't': 0,
    'b': 150,
  }
  return fig

@app.callback(
    Output('G_hist-stats', 'figure'),
    [Input('G_table-stats-signal', 'children'),
     Input('G_table-stats', 'selected_row_indices')])
def update_hist_plot(table_signal, selected_row_indices):
  df = make_table_stats_cache(table_signal)

  # if len(df) <= 5:
    # return ''

  # Determine statistics to plot
  stats_cols = lib.order_chosen_columns(list(df.columns))

  fig = plotly.tools.make_subplots(
    rows = 1, cols = len(stats_cols))

  # Color selected markers
  if len(selected_row_indices) > 0:
    selected_row_index = selected_row_indices[0]
  else:
    selected_row_index = None
  selected_line = dict()

  # Generate each plot
  for idx, stats_col in enumerate(stats_cols):
    subplot_num = idx + 1
    fig.append_trace(
      go.Histogram(
        x = df[stats_col],
        marker = dict(color = lib.get_color(stats_col)),
        name = '',
        opacity = 0.4,
      ), 
      1, subplot_num
    )
    if selected_row_index is not None:
      selected_line[subplot_num] = (df.index[selected_row_index], df[stats_col][len(df) - selected_row_index - 1])

  # Subplot formatting

  all_shapes = []
  x_domains = lib.get_x_domains(len(stats_cols))
  for idx, stats_col in enumerate(stats_cols):
    subplot_num = idx + 1
    fig['layout']['yaxis%s' % (subplot_num)].update(
      fixedrange = True,
      showticklabels = False,
      showgrid = False,
      zeroline = False,
    )
    fig['layout']['xaxis%s' % (subplot_num)].update(
      domain = x_domains[idx],
      title = stats_col,
      fixedrange = True,
      showgrid = True,
      zeroline = False,
      ticks = 'outside',
      ticklen = 3,
      tickcolor = '#eee',
      tickangle = 0, # disable automatic tick rotation
      range = lib.get_batch_statcol_xrange(df[stats_col], stats_col),
    )

    if selected_row_index is not None:
      all_shapes.append(
        lib.get_batch_select_line(
          x0 = selected_line[subplot_num][1],
          x1 = selected_line[subplot_num][1],
          y0 = 0,
          y1 = len(df) / 2.5,
          xref = 'x%s' % (subplot_num),
          yref = 'y1',
        )
      )

  fig['layout']['shapes'] = all_shapes
  # Global figure formatting
  fig['layout']['paper_bgcolor'] = 'rgba(255, 255, 255, 0)'
  fig['layout']['plot_bgcolor'] = 'rgba(255, 255, 255, 0)'
  fig['layout']['showlegend'] = False
  fig['layout']['width'] = 455 + len(stats_cols) * 150
  fig['layout']['height'] = 100
  fig['layout']['margin'] = {
    'l': 430,
    'r': 25,
    't': 0,
    # 't': 60,
    # 'b': 25,
    'b': 40,
  }
  return fig

##
# Download callbacks
##
@app.callback(
  Output('G_download-link', 'href'), 
  [Input('G_hidden-pred-df-stats-signal', 'children')])
def update_link(signal):
  if signal == 'init':
    assert False, 'init'
  stats = grab_s3_stats_cache(signal)

  # Drop extra cols
  drop_cols = [
    '1-bp ins frequency',
    'MH del frequency',
    'MHless del frequency',
  ]
  stats = stats.drop(drop_cols, axis = 1)

  # Rename to shorter versions
  stats = lib.rename_batch_columns(stats)

  # Reformat floats
  stats_cols = list(stats.columns)
  nonstat_cols = [
    'ID',
    'PAM',
    'URL',
    'Cas9 type',
    'Celltype',
    'Chromosome',
    'Cutsite distance to 3p boundary',
    'Cutsite distance to 5p boundary',
    'Exon end',
    # 'Exon number',
    'Exon start',
    'Exon strand',
    'Gene symbol',
    'Genome',
    'Local context',
    'Local cutsite',
    'gRNA',
    'gRNA strand w.r.t. exon strand',
    'kgID',
    'Strand',
    'Cutsite coordinate',
  ]

  for nonstat_col in nonstat_cols:
    stats_cols.remove(nonstat_col)
  for stat_col in stats_cols:
    # Reformat
    if stat_col in ['Precision', 'MH strength']:
      stats[stat_col] = [float('%.2f' % (s)) for s in stats[stat_col]]    
    else:
      stats[stat_col] = [float('%.1f' % (s)) for s in stats[stat_col]]    

  # Reorder columns
  stats = stats[nonstat_cols + lib.order_chosen_columns(stats_cols)]

  time = str(datetime.datetime.now()).replace(' ', '_').replace(':', '-')
  link_fn = '/dash/urlToDownloadGene?value={}'.format(time)
  stats.to_csv('user-csvs/%s.csv' % (time), index = False)
  return link_fn

@app.callback(
  Output('G_download-link', 'children'), 
  [Input('G_hidden-pred-df-stats-signal', 'children')])
def update_link_text(signal):
  if signal == 'init':
    assert False, 'init'
  stats = grab_s3_stats_cache(signal)
  num_grnas = len(stats)
  num_kgids = len(set(stats['kgID']))
  return '📑 Download full table of predictions for %s gRNAs and %s kgIDs' % (num_grnas, num_kgids)

##
# Flask serving
##
@app.server.route('/dash/urlToDownloadGene') 
def download_csv_gene():
  value = flask.request.args.get('value')
  # create a dynamic csv or file here using `StringIO` 
  # (instead of writing to the file system)
  local_csv_fn = value.split('/')[-1]
  return flask.send_file(
    open('user-csvs/%s.csv' % (local_csv_fn), 'rb'),
    mimetype = 'text/csv',
    attachment_filename = 'inDelphi_gene_output.csv',
    as_attachment = True,
  )


##
# Page link callback
##
@app.callback(
  Output('G_page-link', 'href'),
  [Input('G_genome-radio', 'value'),
   Input('G_gene-dropdown', 'value'),
   Input('G_celltype_dropdown', 'value'),
   Input('G_dropdown-columns', 'value'),
   Input('G_dropdown-columns', 'options'),
   Input('G_dropdown-sortcol', 'value'),
   Input('G_sortdirection', 'value'),
   Input('G_table-stats', 'selected_row_indices'),
  ])
def update_pagelink(genome_build, gene, celltype, chosen_columns, column_options, sort_by, sort_dir, selected_row):
  url = '%s' % (lib.encode_url_path_gene(genome_build, gene, celltype, chosen_columns, column_options, sort_by, sort_dir, selected_row))
  return url