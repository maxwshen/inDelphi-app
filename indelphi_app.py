import pickle, copy, os, datetime, subprocess
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

import inDelphi
import generalStats
import lib

# Dash server setup
app = dash.Dash('')
server = app.server

# init
inDelphi.init_model()
if not os.path.isdir('user-csvs/'):
  os.mkdir('user-csvs/')
else:
  subprocess.check_output('rm -rf user-csvs/*', shell = True)


# Remove these plotly modebar buttons to limit interactivity
modebarbuttons_2d = ['zoom2d', 'pan2d', 'select2d', 'lasso2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'resetScale2d', 'hoverClosestCartesian', 'hoverCompareCartesian', 'toggleSpikelines']

##
headerHeight = 130


###################################################################
###################################################################
##
# App layout
##
app.layout = html.Div([
  html.Div([

    ##
    # Hidden divs for light data storage
    ##
    html.Div(
      [
        html.Div(
          id = 'hidden-pred-df',
          children = 'init'
        ),
        html.Div(
          id = 'hidden-pred-stats',
          children = 'init'
        ),
        html.Div(
          id = 'hidden-cache-dsb-left',
          children = 'init'
        ),
        html.Div(
          id = 'hidden-cache-dsb-right',
          children = 'init'
        ),
        html.Div(
          id = 'hidden-cache-pam-left',
          children = 'init'
        ),
        html.Div(
          id = 'hidden-cache-pam-right',
          children = 'init'
        ),

        dcc.Location(
          id = 'url',
          refresh = False,
        ),
      ],
      style = dict(
        display = 'none',
      ),
    ),

    ##
    # Header
    ##
    html.Div(
      [
        ###################################################
        # Upper header
        ###################################################
        html.H4(
          'inDelphi',
          style = dict(
            textAlign = 'center',
          ),
        ),

        ###################################################
        # Sequence boxes
        ###################################################
        html.Div([
          # Left box
          html.Div(
            [
              dcc.Input(
                id = 'textbox1', 
                size = 28,
                value = 'TAGTTTCTAGCACAGCGCTGGTGTGG',
                type = 'text',
                autofocus = True,
                style = dict(
                  textAlign = 'right',
                  direction = 'rtl',
                  fontFamily = 'monospace',
                  fontSize = 16,
                  float = 'right',
                ),
              )
            ],
            className = 'dna_textbox',
          ),

          # Right box
          html.Div(
            [
              dcc.Input(
                id = 'textbox2', 
                size = 28,
                value = 'CGTGTGGCTGAAGGCATAGTAATTCTGA',
                type = 'text',
                style = dict(
                  textAlign = 'left',
                  fontFamily = 'monospace',
                  fontSize = 16,
                  float = 'left',
                ),
              ),
            ],
            className = 'dna_textbox',
          ),
        ], 
          style = dict(
            verticalAlign = 'center',
            whiteSpace = 'nowrap',
            overflowX = 'auto',
          ),
        ),

        ###################################################
        # Row of PAM arrows 
        ###################################################
        html.Div(
          [
            html.Div(
              [
                # (Left) placeholder
                html.Div('',
                  style = dict(
                    display = 'table-cell',
                    width = '25%',
                  ),
                ),

                ###################################################
                # (Center) DSB arrows: Simple
                ###################################################
                html.Div([
                    html.A('◄',
                      id = 'button-dsb-left',
                      style = dict(
                        textDecoration = 'none',
                        fontFamily = 'monospace',
                        color = 'rgb(0, 174, 179)',
                        fontSize = 20,
                      ),
                    ),
                    '\tDSB\t',
                    html.A('►',
                      id = 'button-dsb-right',
                      style = dict(
                        textDecoration = 'none',
                        fontFamily = 'monospace',
                        color = 'rgb(0, 174, 179)',
                        fontSize = 20,
                      ),
                    ),
                  ],
                  style = dict(
                    textAlign = 'center',
                    display = 'table-cell',
                    width = '50%',
                  )
                ),

                ###################################################
                # (Right) PAM arrows
                ###################################################
                html.Div([
                    html.A('◄',
                      id = 'button-pam-left',
                      style = dict(
                        textDecoration = 'none',
                        fontFamily = 'monospace',
                        fontSize = 20,
                        color = 'rgb(237, 71, 149)',
                        verticalAlign = 'middle',
                      ),
                    ),
                    dcc.Input(
                      id = 'textbox_pam', 
                      size = 5,
                      value = 'NGG',
                      type = 'text',
                      autofocus = True,
                      style = dict(
                        fontFamily = 'monospace',
                        fontSize = 14,
                        textAlign = 'center',
                        height = '18px',
                        width = '70px',
                        marginLeft = '5px',
                        marginRight = '5px',
                      ),
                    ),
                    html.A('►',
                      id = 'button-pam-right',
                      style = dict(
                        textDecoration = 'none',
                        fontFamily = 'monospace',
                        fontSize = 20,
                        color = 'rgb(237, 71, 149)',
                        verticalAlign = 'middle',
                      ),
                    ),
                  ],
                  style = dict(
                    display = 'table-cell',
                    textAlign = 'right',
                    width = '25%',
                  )
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
          ),
        ),



      ],
      style = dict(
        position = 'fixed',
        backgroundColor = 'white',
        borderBottom = '3px solid #777777',
        zIndex = 1e6,
        width = '1010px',
        left = '50%',
        transform = 'translate(-50%, 0)',
        height = headerHeight,
        marginTop = '-%spx' % (headerHeight + 20),
      ),
    ),

    ##
    # Body / plots
    ##
    html.Div(
      [
        ###################################################
        # Module: Summary, alignment preview
        ###################################################
        html.Div([
          # header
          html.Div([
            html.Div([
              html.Strong('Summary of predictions: Top 10 frequent events')
              ],
              className = 'module_header_text'),
            ],
            className = 'module_header'
          ),

          html.Div(
            [
              # Text table
              dcc.Graph(
                id = 'summary-alignment-table',
                config = dict(
                  modeBarButtonsToRemove = modebarbuttons_2d,
                  displaylogo = False,
                ),
                style = dict(
                  height = 300,
                ),
                className = 'eight columns',
              ),

              # bar chart
              dcc.Graph(
                id = 'summary-alignment-barchart',
                config = dict(
                  modeBarButtonsToRemove = modebarbuttons_2d,
                  displaylogo = False,
                ),
                style = dict(
                  height = 300,
                ),
                className = 'four columns',
              ),
            ],
            className = 'row',
          ),
        ], className = 'module_style',
        ),


        ###################################################
        # Module: Indel Lengths
        ###################################################
        html.Div([
          # header
          html.Div([
            html.Div([
              html.Strong('Indel length predictions')
              ],
              className = 'module_header_text'),
            ],
            className = 'module_header'
          ),

          html.Div(
            [
              # Frameshift
              html.Div(
                [
                  dcc.Graph(
                    id = 'plot-fs',
                    style = dict(
                      height = 300, 
                    ),
                    config = dict(
                      modeBarButtonsToRemove = modebarbuttons_2d,
                      displaylogo = False,
                    ),
                  ),
                ],
                className = 'three columns',
              ),

              # Indel length
              html.Div(
                [
                  dcc.Graph(
                    id = 'plot-indel-len',
                    style = dict(
                      height = 300, 
                    ),
                    config = dict(
                      modeBarButtonsToRemove = modebarbuttons_2d,
                      displaylogo = False,
                    ),
                  ),
                ],
                className = 'nine columns',
              ),
            ],
            className = 'row',
          ),
        ], className = 'module_style',
        ),

        ###################################################
        # Module: Genome statistics
        ###################################################
        ## Precision
        html.Div([
          # header
          html.Div([
            html.Div([
              html.Strong('Comparison to predictions at {:,} SpCas9 target sites in human exons and introns'.format(int(sum(generalStats.gs_logphi.unnormY))))
              ],
              className = 'module_header_text'),
            ],
            className = 'module_header'
          ),

          html.Div(
            [
              html.Div(
                [
                  dcc.Graph(
                    id = 'plot-genstats-precision',
                    style = dict(
                      height = 200, 
                      width = 970/2,
                    ),
                    config = dict(
                      modeBarButtonsToRemove = modebarbuttons_2d,
                      displaylogo = False,
                    ),
                  ),
                ],
                className = 'six columns',
              ),
              html.Div(
                [
                  html.Div(
                    id = 'text-genstats-precision',
                    className = 'generalstats_text_inner',
                  ),
                ],
                style = dict(
                  height = 200,
                ),
                className = 'six columns',
              ),
            ],
            className = 'row',
          ),

          ## Phi
          html.Div(
            [
              html.Div(
                [
                  dcc.Graph(
                    id = 'plot-genstats-logphi',
                    style = dict(
                      height = 200, 
                      width = 970/2,
                    ),
                    config = dict(
                      modeBarButtonsToRemove = modebarbuttons_2d,
                      displaylogo = False,
                    ),
                  ),
                ],
                className = 'six columns',
              ),
              html.Div(
                [
                  html.Div(
                    id = 'text-genstats-logphi',
                    className = 'generalstats_text_inner',
                  ),
                ],
                style = dict(
                  height = 200,
                ),
                className = 'six columns',
              ),
            ],
            className = 'row',
          ),

          ## Frameshift frequency
          html.Div(
            [
              html.Div(
                [
                  dcc.Graph(
                    id = 'plot-genstats-frameshift',
                    style = dict(
                      height = 200, 
                      width = 970/2,
                    ),
                    config = dict(
                      modeBarButtonsToRemove = modebarbuttons_2d,
                      displaylogo = False,
                    ),
                  ),
                ],
                className = 'six columns',
              ),
              html.Div(
                [
                  html.Div(
                    id = 'text-genstats-frameshift',
                    className = 'generalstats_text_inner',
                  ),
                ],
                style = dict(
                  height = 200,
                ),
                className = 'six columns',
              ),
            ],
            className = 'row',
          ),
        ], className = 'module_style',
        ),

        ###################################################
        # Module: Detailed genotypes
        ###################################################
        html.Div([
          # header
          html.Div([
            html.Div([
              html.Strong('All predictions of 1-bp insertion and 1- to 60-bp deletion events')
              ],
              className = 'module_header_text'),
            ],
            className = 'module_header'
          ),

          dcc.Graph(
            id = 'plot-table-genotypes',
            style = dict(
            ),
            config = dict(
              modeBarButtonsToRemove = modebarbuttons_2d,
              displaylogo = False,
            ),
          ),

          dt.DataTable(
            id = 'table-genotypes',
            rows = [{}], # init rows
            row_selectable = True,
            filterable = True,
            sortable = True,
            selected_row_indices = [],
          ),

          html.Div([
            html.Div(
              html.A(
                'Download CSV of inDelphi predictions', 
                id = 'csv-download-link'
              ),
            ),
            html.Div(
              html.A(
                'Shareable link to your results', 
                id = 'page-link'
              ),
            ),
            ],
            style = dict(
              marginLeft = '20',
            ),
          ),

          html.Div(
            'Copyright MIT 2018.\nAll Rights Reserved.',
            style = dict(
              textAlign = 'center',
              marginTop = '30',
              marginBottom = '30',
            )
          ),
        ], className = 'module_style',
        ),

      ],
      # body style
      style = dict(
        marginTop = '%spx' % (headerHeight + 20),
      ),
    ),
    ##
  ], 
    style = dict(
      width = '970px',
      margin = '0 auto',
    )
  ),
],  # body div
style = dict(
  width = '1000px',
  margin = '0 auto',
)
)

#######################################################################
#########################      CALLBACKS      #########################
#######################################################################
##
# Header Callbacks
##

## Arrow buttons
# These must be 1->1 mapping, otherwise we can't tell which n_clicks triggered the callback
@app.callback(
  Output('hidden-cache-dsb-left', 'children'),
  [Input('button-dsb-left', 'n_clicks')])
def cb_update_cache_dsb_left(n_clicks):
  return '%s' % (time.time())

@app.callback(
  Output('hidden-cache-dsb-right', 'children'),
  [Input('button-dsb-right', 'n_clicks')])
def cb_update_cache_dsb_right(n_clicks):
  return '%s' % (time.time())

@app.callback(
  Output('hidden-cache-pam-left', 'children'),
  [Input('button-pam-left', 'n_clicks')])
def cb_update_cache_pam_left(n_clicks):
  return '%s' % (time.time())

@app.callback(
  Output('hidden-cache-pam-right', 'children'),
  [Input('button-pam-right', 'n_clicks')])
def cb_update_cache_pam_right(n_clicks):
  return '%s' % (time.time())


@app.callback(
  Output('textbox1', 'value'),
  [Input('hidden-cache-dsb-left', 'children'),
   Input('hidden-cache-dsb-right', 'children'),
   Input('hidden-cache-pam-left', 'children'),
   Input('hidden-cache-pam-right', 'children'),
   Input('url', 'pathname')],
  [State('textbox1', 'value'),
   State('textbox2', 'value'),
   State('textbox_pam', 'value')])
def cb_update_textbox1_arrow(cache_dsb_left, cache_dsb_right, cache_pam_left, cache_pam_right, url, text1, text2, text_pam):
  left_dsb_time = float(cache_dsb_left)
  right_dsb_time = float(cache_dsb_right)
  pageload_dsb = bool(abs(left_dsb_time - right_dsb_time) < 0.01)
  left_pam_time = float(cache_pam_left)
  right_pam_time = float(cache_pam_right)
  pageload_pam = bool(abs(left_pam_time - right_pam_time) < 0.01)
  if pageload_dsb and pageload_pam:
    valid_flag, seq, cutsite = lib.parse_valid_url_path(url)
    if not valid_flag or cutsite is None:
      return text1
    else:
      return seq[:cutsite]

  if max(left_dsb_time, right_dsb_time) > max(left_pam_time, right_pam_time):
    # Clicked DSB
    if left_dsb_time > right_dsb_time:
      return text1[:-1]
    elif right_dsb_time > left_dsb_time:
      return text1 + text2[0]
  elif max(left_pam_time, right_pam_time) > max(left_dsb_time, right_dsb_time):
    # Clicked PAM
    if left_pam_time > right_pam_time:
      ## TO IMPLEMENT
      # Seek left PAM using text_pam
      return text1[:-1]
    elif right_pam_time > left_pam_time:
      ## TO IMPLEMENT
      # Seek right PAM using text_pam
      return text1 + text2[0]

@app.callback(
  Output('textbox2', 'value'),
  [Input('hidden-cache-dsb-left', 'children'),
   Input('hidden-cache-dsb-right', 'children'),
   Input('hidden-cache-pam-left', 'children'),
   Input('hidden-cache-pam-right', 'children'),
   Input('url', 'pathname')],
  [State('textbox1', 'value'),
   State('textbox2', 'value'),
   State('textbox_pam', 'value')])
def cb_update_textbox2_arrow(cache_dsb_left, cache_dsb_right, cache_pam_left, cache_pam_right, url, text1, text2, text_pam):
  left_dsb_time = float(cache_dsb_left)
  right_dsb_time = float(cache_dsb_right)
  pageload_dsb = bool(abs(left_dsb_time - right_dsb_time) < 0.01)
  left_pam_time = float(cache_pam_left)
  right_pam_time = float(cache_pam_right)
  pageload_pam = bool(abs(left_pam_time - right_pam_time) < 0.01)
  if pageload_dsb and pageload_pam:
    valid_flag, seq, cutsite = lib.parse_valid_url_path(url)
    if not valid_flag or cutsite is None:
      return text2
    else:
      return seq[cutsite:]

  if max(left_dsb_time, right_dsb_time) > max(left_pam_time, right_pam_time):
    # Clicked DSB
    if left_dsb_time > right_dsb_time:
      return text1[-1] + text2
    elif right_dsb_time > left_dsb_time:
      return text2[1:]
  elif max(left_pam_time, right_pam_time) > max(left_dsb_time, right_dsb_time):
    # Clicked PAM
    if left_dsb_time > right_dsb_time:
      ## TO IMPLEMENT
      return text1[-1] + text2
    elif right_dsb_time > left_dsb_time:
      ## TO IMPLEMENT
      return text2[1:]


##
# Prediction callback
##
@app.callback(
  Output('hidden-pred-df', 'children'),
  [Input('textbox1', 'value'),
   Input('textbox2', 'value')])
def cb_update_pred_df(text1, text2):
  seq = text1 + text2
  cutsite = len(text1)
  pred_df, stats = inDelphi.predict(seq, cutsite)
  return pred_df.to_csv()

@app.callback(
  Output('hidden-pred-stats', 'children'),
  [Input('textbox1', 'value'),
   Input('textbox2', 'value')])
def cb_update_pred_stats(text1, text2):
  seq = text1 + text2
  cutsite = len(text1)
  pred_df, stats = inDelphi.predict(seq, cutsite)
  return pd.DataFrame(stats, index = [0]).to_csv()

##
# Summary of predictions callbacks
##
@app.callback(
  Output('summary-alignment-table', 'figure'),
  [Input('hidden-pred-df', 'children'),
   Input('hidden-pred-stats', 'children'),
  ])
def cb_update_summary_alignment_text(pred_df_string, pred_stats_string):
  pred_df = pd.read_csv(StringIO(pred_df_string), index_col = 0)
  stats = pd.read_csv(StringIO(pred_stats_string), index_col = 0)
  
  inDelphi.add_genotype_column(pred_df, stats)

  top10 = pred_df.sort_values('Predicted frequency', ascending = False).iloc[:10]
  gts = top10['Genotype']
  fqs = top10['Predicted frequency']
  lens = top10['Length']
  cats = top10['Category']
  fq_strings = ['-'] + ['%.1f' % (s) for s in fqs]

  def trim_alignment(gt, cutsite, name):
    radius = 26
    if name == 'ins':
      trim_cand = gt[cutsite - radius : cutsite + radius + 1]
      if len(trim_cand) == 2*radius + 1:
        return trim_cand
      else:
        return gt
    else:
      trim_cand = gt[cutsite - radius : cutsite + radius + 1]
      if len(trim_cand) == 2*radius + 1:
        return trim_cand
      else:
        return gt
    return

  def add_bar(seq, cutsite):
    return seq[:cutsite] + '|' + seq[cutsite:]

  def get_gapped_alignments(top, stats):
    cutsite = stats['Cutsite'].iloc[0]
    gapped_aligns = []
    for idx, row in top.iterrows():
      gt = row['Genotype']
      gt_pos = row['Genotype position']
      length = row['Length']
      cat = row['Category']
      if cat == 'ins':
        gapped_aligns.append(trim_alignment(gt, cutsite, 'ins'))
        continue
      if gt_pos == 'e':
        gapped_aligns.append('multiple deletion genotypes')
        continue

      gt_pos = int(gt_pos)
      gap_gt = gt[:cutsite - length + gt_pos] + '-'*length + gt[cutsite - length + gt_pos:]
      gap_gt = add_bar(gap_gt, cutsite)
      gapped_aligns.append(trim_alignment(gap_gt, cutsite, 'del'))
    return gapped_aligns

  gap_gts = get_gapped_alignments(top10, stats)
  
  alignments, categories = [], []
  cutsite = stats['Cutsite'].iloc[0]
  reference_seq = stats['Reference sequence'].iloc[0]
  reference_seq = add_bar(reference_seq, cutsite)
  alignments.append(trim_alignment(reference_seq, cutsite, 'ref'))
  categories.append('Reference')

  for gt, length, cat in zip(gap_gts, lens, cats):
    alignments.append(gt)
    if cat == 'ins':
      categories.append('%s-bp insertion' % (length))
    elif cat == 'del':
      categories.append('%s-bp deletion' % (length))

  return dict(
    data = [go.Table(
      type = 'table',
      columnwidth = [300, 100, 80],
      header = dict(
        values = ['Alignment', 'Category', '%'],
        align = ['center', 'right', 'right'], 
        line = dict(color = 'white'),
        fill = dict(color = 'white'),
      ),
      cells = dict(
        values = [alignments, categories, fq_strings],
        align = ['center', 'right', 'right'], 
        line = dict(color = 'white'),
        fill = dict(color = 'white'),
        font = dict(family = 'monospace'),
      )
    )],
    layout = go.Layout(
      font = dict(
        family = 'monospace',
      ),
      margin = dict(
        l = 10,
        r = 0,
        t = 5,
        b = 5,
      ),
    ),
  )

@app.callback(
  Output('summary-alignment-barchart', 'figure'),
  [Input('hidden-pred-df', 'children'),
   Input('hidden-pred-stats', 'children'),
  ])
def cb_update_summary_alignment_barchart(pred_df_string, pred_stats_string):
  pred_df = pd.read_csv(StringIO(pred_df_string), index_col = 0)
  stats = pd.read_csv(StringIO(pred_stats_string), index_col = 0)
  
  inDelphi.add_genotype_column(pred_df, stats)

  top10 = pred_df.sort_values('Predicted frequency', ascending = False).iloc[:10]
  fqs = top10['Predicted frequency'][::-1]

  cats = top10['Category'][::-1]
  gts = top10['Genotype position'][::-1]
  colors = []
  for cat, gt in zip(cats, gts):
    if gt == 'e':
      colors.append('rgb(236, 100, 12)')
    else:
      if cat == 'del':
        colors.append('rgb(221, 46, 31)')
      else:
        colors.append('rgb(0, 160, 220)')

  return dict(
    data = [go.Bar(
      x = fqs,
      y = np.arange(len(fqs)),
      orientation = 'h',
      marker = dict(
        color = colors,
        line = dict(
          width = 0,
        ),
      ),
    )],
    layout = go.Layout(
      yaxis = dict(
        showticklabels = False,
      ),
      margin = dict(
        l = 5,
        r = 10,
        t = 58,
        b = 42,
      ),
    ),
  )

##
# General stats callbacks
##
@app.callback(
  Output('plot-genstats-precision', 'figure'),
  [Input('hidden-pred-df', 'children'),
   Input('hidden-pred-stats', 'children'),
  ])
def cb_plot_genstats_precision(pred_df_string, pred_stats_string):
  pred_df = pd.read_csv(StringIO(pred_df_string), index_col = 0)
  stats = pd.read_csv(StringIO(pred_stats_string), index_col = 0)
  xval = stats['Precision'].iloc[0]
  return dict(
    data = [
      generalStats.gs_precision.trace(xval),
    ],
    layout = generalStats.gs_precision.layout(xval),
  )

@app.callback(
  Output('plot-genstats-logphi', 'figure'),
  [Input('hidden-pred-df', 'children'),
   Input('hidden-pred-stats', 'children'),
  ])
def cb_plot_genstats_logphi(pred_df_string, pred_stats_string):
  pred_df = pd.read_csv(StringIO(pred_df_string), index_col = 0)
  stats = pd.read_csv(StringIO(pred_stats_string), index_col = 0)
  xval = np.log(stats['Phi'].iloc[0])
  return dict(
    data = [
      generalStats.gs_logphi.trace(xval),
    ],
    layout = generalStats.gs_logphi.layout(xval),
  )

@app.callback(
  Output('plot-genstats-frameshift', 'figure'),
  [Input('hidden-pred-df', 'children'),
   Input('hidden-pred-stats', 'children'),
  ])
def cb_plot_genstats_frameshift(pred_df_string, pred_stats_string):
  pred_df = pd.read_csv(StringIO(pred_df_string), index_col = 0)
  stats = pd.read_csv(StringIO(pred_stats_string), index_col = 0)
  xval = stats['Frameshift frequency'].iloc[0]
  return dict(
    data = [
      generalStats.gs_frameshift.trace(xval),
    ],
    layout = generalStats.gs_frameshift.layout(xval),
  )


## General stats text
@app.callback(
  Output('text-genstats-precision', 'children'),
  [Input('hidden-pred-df', 'children'),
   Input('hidden-pred-stats', 'children'),
  ])
def cb_text_genstats_precision(pred_df_string, pred_stats_string):
  pred_df = pd.read_csv(StringIO(pred_df_string), index_col = 0)
  stats = pd.read_csv(StringIO(pred_stats_string), index_col = 0)
  xval = stats['Precision'].iloc[0]
  cum, var_text, var_color = generalStats.gs_precision.cumulative(xval)
  return [
    html.Strong('This target site has '),
    html.Strong(var_text,
      style = dict(color = var_color),
    ),
    html.Strong(' precision.',
      style = dict(color = var_color),
    ),
    html.Br(),
    html.Span('Precision score: %.2f' % (xval),
      className = 'generalstats_subtext_style'),
    html.Br(),
    html.Span('Percentile: %s' % (cum),
      className = 'generalstats_subtext_style'),
  ]

@app.callback(
  Output('text-genstats-logphi', 'children'),
  [Input('hidden-pred-df', 'children'),
   Input('hidden-pred-stats', 'children'),
  ])
def cb_text_genstats_logphi(pred_df_string, pred_stats_string):
  pred_df = pd.read_csv(StringIO(pred_df_string), index_col = 0)
  stats = pd.read_csv(StringIO(pred_stats_string), index_col = 0)
  xval = np.log(stats['Phi'].iloc[0])
  cum, var_text, var_color = generalStats.gs_logphi.cumulative(xval)
  return [
    html.Strong('This target site has '),
    html.Strong(var_text,
      style = dict(color = var_color),
    ),
    html.Strong(' microhomology strength.',
      style = dict(color = var_color),
    ),
    html.Br(),
    html.Span('Log phi: %.2f' % (xval),
      className = 'generalstats_subtext_style'),
    html.Br(),
    html.Span('Percentile: %s' % (cum),
      className = 'generalstats_subtext_style'),
  ]

@app.callback(
  Output('text-genstats-frameshift', 'children'),
  [Input('hidden-pred-df', 'children'),
   Input('hidden-pred-stats', 'children'),
  ])
def cb_text_genstats_frameshift(pred_df_string, pred_stats_string):
  pred_df = pd.read_csv(StringIO(pred_df_string), index_col = 0)
  stats = pd.read_csv(StringIO(pred_stats_string), index_col = 0)
  xval = stats['Frameshift frequency'].iloc[0]
  cum, var_text, var_color = generalStats.gs_frameshift.cumulative(xval)
  return [
    html.Strong('This target site has '),
    html.Strong(var_text,
      style = dict(color = var_color),
    ),
    html.Strong(' frameshift frequency.',
      style = dict(color = var_color),
    ),
    html.Br(),
    html.Span('Frameshift frequency: %.1f%%' % (xval),
      className = 'generalstats_subtext_style'),
    html.Br(),
    html.Span('Percentile: %s' % (cum),
      className = 'generalstats_subtext_style'),
  ]

##
# Indel length and frameshift callbacks
@app.callback(
  Output('plot-indel-len', 'figure'),
  [Input('hidden-pred-df', 'children'),
  ])
def cb_plot_indel_len(pred_df_string):
  pred_df = pd.read_csv(StringIO(pred_df_string), index_col = 0)

  lendf = inDelphi.get_indel_length_fqs(pred_df)

  X = [int(s) for s in lendf['Indel length']]
  Y = [s for s in lendf['Predicted frequency']]
  return dict(
    data = [
      go.Bar(
        x = X,
        y = Y,
        opacity = 0.6,
        width = 0.9, 
        marker = dict(
          color = 'rgb(200, 20, 20)',
          line = dict(
            width = 0,
          ),
        ),
      )
    ],
    layout = go.Layout(
      xaxis = dict(
        autorange = 'reversed',
        title = 'Indel length',
        ticks = 'outside',
        ticklen = 3,
        tickwidth = 0.5,
        tick0 = 0,
        dtick = 5,
        zeroline = False,
      ),
      yaxis = dict(
        title = 'Frequency (%)',
        hoverformat = '.2f%%',
        zeroline = False,
      ),
      font = dict(
        family = 'Arial',
      ),
      margin = dict(
        t = 30,
        l = 70,
        r = 20,
        b = 70,
      ),
    ),
  )

@app.callback(
  Output('plot-fs', 'figure'),
  [Input('hidden-pred-df', 'children'),
  ])
def cb_plot_fs(pred_df_string):
  pred_df = pd.read_csv(StringIO(pred_df_string), index_col = 0)

  fs_df = inDelphi.get_frameshift_fqs(pred_df)
  X = ['+0', '+1', '+2']
  Y = [float(fs_df[fs_df['Frame'] == s]['Predicted frequency']) for s in X]
  return dict(
    data = [
      go.Bar(
        x = X,
        y = Y,
        text = ['%.0f%%' % (s) for s in Y],
        textfont = dict(
          family = 'Arial',
        ),
        textposition = 'auto',
        opacity = 0.6,
        # width = 1,  # touching
        marker = dict(
          color = 'rgb(158, 202, 225)',
          line = dict(
            color = 'rgb(8, 48, 107)',
            width = 1.5,
          ),
        ),
      ),
    ],
    layout = go.Layout(
      # title = 'Frameshift frequency (%)',
      # margin = go.Margin(l = 40, r = 0, t = 40, b = 30),
      xaxis = dict(
        title = 'Frame',
        showline = False,
        tickvals = list(range(3)),
        ticktext = ['+0', '+1', '+2'],
      ),
      yaxis = dict(
        title = 'Frameshift frequency (%)',
        titlefont = dict(
          family = 'Arial',
        ),
        zeroline = False,
        showline = True,
        ticks = 'outside',
        tick0 = 0,
        ticklen = 3,
        tickwidth = 0.5,
        hoverformat = '.2f%%',
      ),
      font = dict(
        family = 'Arial',
      ),
      margin = dict(
        t = 30,
        l = 70,
        r = 20,
        b = 70,
      ),
    ),
  )



##
# Genotype table callbacks
## 
@app.callback(
  Output('table-genotypes', 'rows'), 
  [Input('hidden-pred-df', 'children'),
   Input('hidden-pred-stats', 'children'),
  ])
def cb_update_genotype_table(pred_df_string, pred_stats_string):
  pred_df = pd.read_csv(StringIO(pred_df_string), index_col = 0)
  stats = pd.read_csv(StringIO(pred_stats_string), index_col = 0)

  inDelphi.add_genotype_column(pred_df, stats)
  inDelphi.add_name_column(pred_df, stats)

  # Edit for display
  pred_df['Frequency (%)'] = ['%.1f' % (s) for s in pred_df['Predicted frequency']]
  pred_df = pred_df.drop(
    [
      'Inserted Bases', 
      'Genotype position',
      'Predicted frequency',
      'Category',
    ], 
    axis = 1
  )
  pred_df = pred_df[['Name', 'Length', 'Frequency (%)', 'Genotype']]
  return pred_df.to_dict('records')

@app.callback(
  Output('plot-table-genotypes', 'figure'),
  [Input('table-genotypes', 'rows'),
   Input('table-genotypes', 'selected_row_indices')
  ])
def cb_update_genotype_plot(rows, selected_row_indices):
  df = pd.DataFrame(rows)
  colors =  ['#0074D9'] * len(df)
  for idx in (selected_row_indices or []):
    colors[idx] = '#FF851B'
  return dict(
    data = [
      go.Bar(
        x = df['Name'],
        y = df['Frequency (%)'],
        marker = dict(
          color = colors,
          line = dict(
            width = 0,
          ),
        ),
      ),
    ],
    layout = go.Layout(
      xaxis = dict(
      ),
      yaxis = dict(
        title = 'Frequency (%)',
      ),
      font = dict(
        family = 'Arial',
      ),
      margin = dict(
        t = 10,
        r = 20,
      ),
    ),
  )

@app.callback(
  Output('table-genotypes', 'selected_row_indices'),
  [Input('plot-table-genotypes', 'clickData')],
  [State('table-genotypes', 'selected_row_indices')]
  )
def cb_update_datatable_selected(clickData, selected_row_indices):
  # Update selections in table based on clicking plot
  if clickData:
    for point in clickData['points']:
      if point['pointNumber'] in selected_row_indices:
        selected_row_indices.remove(point['pointNumber'])
      else:
        selected_row_indices.append(point['pointNumber'])
  return selected_row_indices


##
# Download callbacks
##
@app.callback(
  Output('csv-download-link', 'href'), 
  [Input('hidden-pred-df', 'children'),
   Input('hidden-pred-stats', 'children'),
  ])
def cb_update_link(pred_df_string, pred_stats_string):
  pred_df = pd.read_csv(StringIO(pred_df_string), index_col = 0)
  stats = pd.read_csv(StringIO(pred_stats_string), index_col = 0)

  inDelphi.add_genotype_column(pred_df, stats)
  inDelphi.add_name_column(pred_df, stats)

  time = str(datetime.datetime.now()).replace(' ', '_').replace(':', '-')
  link_fn = '/dash/urlToDownload?value={}'.format(time)
  pred_df.to_csv('user-csvs/%s.csv' % (time))

  return link_fn

@app.server.route('/dash/urlToDownload') 
def download_csv():
  value = flask.request.args.get('value')
  # create a dynamic csv or file here using `StringIO` 
  # (instead of writing to the file system)
  local_csv_fn = value.split('/')[-1]
  return flask.send_file(
    open('user-csvs/%s' % (local_csv_fn + '.csv'), 'rb'),
    mimetype = 'text/csv',
    attachment_filename = 'inDelphi_output.csv',
    as_attachment = True,
  )

##
# Page link callback
##
@app.callback(
  Output('page-link', 'href'),
  [Input('textbox1', 'value',),
   Input('textbox2', 'value',),
  ])
def cb_update_pagelink(text1, text2):
  seq = text1 + text2
  cutsite = len(text1)
  return 'https://dev.crisprindelphi.design/%s' % (lib.encode_dna_to_url_path(seq, cutsite))

##
# Local CSS
##
css_directory = os.getcwd()
stylesheets = ['stylesheet.css']
@app.server.route('/static/<stylesheet>')
def serve_stylesheet(stylesheet):
  if stylesheet not in stylesheets:
    raise Exception(
      '"{}" is excluded from the allowed static files'.format(
        stylesheet
      )
    )
  return flask.send_from_directory(css_directory, stylesheet)


###################################################################
###################################################################
# CSS
for stylesheet in stylesheets:
  app.css.append_css({'external_url': '/static/{}'.format(stylesheet)})

if __name__ == '__main__':
  app.run_server()