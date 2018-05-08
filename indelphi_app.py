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
          id = 'hidden-cache-left',
          children = 'init'
        ),
        html.Div(
          id = 'hidden-cache-right',
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
        html.H4(
          'inDelphi',
          style = dict(
            textAlign = 'center',
          ),
        ),

        # Row
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

        html.Div([
            html.A('◄',
              id = 'button-input-left',
              style = dict(
                textDecoration = 'none',
                fontFamily = 'monospace',
                fontSize = 20,
              ),
            ),
            '\tDSB\t',
            html.A('►',
              id = 'button-input-right',
              style = dict(
                textDecoration = 'none',
                fontFamily = 'monospace',
                fontSize = 20,
              ),
            ),
          ],
          style = dict(
            textAlign = 'center',
          )
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
              html.Strong('Summary of predictions')
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
                    id = 'summary-plot-piechart',
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
                id = 'summary-alignment-text',
                className = 'nine columns',
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
@app.callback(
  Output('hidden-cache-left', 'children'),
  [Input('button-input-left', 'n_clicks')],
  [State('textbox1', 'value')])
def cb_update_cache_left(n_clicks, box1_val):
  return '%s_%s' % (box1_val[-1], time.time())

@app.callback(
  Output('hidden-cache-right', 'children'),
  [Input('button-input-right', 'n_clicks')],
  [State('textbox2', 'value')])
def cb_update_cache_right(n_clicks, box2_val):
  return '%s_%s' % (box2_val[0], time.time())

@app.callback(
  Output('textbox1', 'value'),
  [Input('hidden-cache-left', 'children'),
   Input('hidden-cache-right', 'children'),
   Input('url', 'pathname')],
  [State('textbox1', 'value')])
def cb_update_textbox1_arrow(cache_left, cache_right, url, text):
  left_char = cache_left.split('_')[0]
  left_time = float(cache_left.split('_')[1])
  right_char = cache_right.split('_')[0]
  right_time = float(cache_right.split('_')[1])
  if abs(left_time - right_time) < 0.01:
    valid_flag, seq, cutsite = lib.parse_valid_url_path(url)
    if not valid_flag or cutsite is None:
      return text
    else:
      return seq[:cutsite]
  if left_time > right_time:
    return text[:-1]
  else:
    return text + right_char

@app.callback(
  Output('textbox2', 'value'),
  [Input('hidden-cache-left', 'children'),
   Input('hidden-cache-right', 'children'),
   Input('url', 'pathname')],
  [State('textbox2', 'value')])
def cb_update_textbox2_arrow(cache_left, cache_right, url, text):
  left_char = cache_left.split('_')[0]
  left_time = float(cache_left.split('_')[1])
  right_char = cache_right.split('_')[0]
  right_time = float(cache_right.split('_')[1])
  if abs(left_time - right_time) < 0.01:
    valid_flag, seq, cutsite = lib.parse_valid_url_path(url)
    if not valid_flag or cutsite is None:
      return text
    else:
      return seq[cutsite:]
  if left_time > right_time:
    return left_char + text
  else:
    return text[1:]

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
  Output('summary-alignment-text', 'children'),
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

  def get_gapped_alignments(top, stats):
    cutsite = stats['Cutsite'].iloc[0]
    gapped_aligns = []
    for idx, row in top.iterrows():
      gt = row['Genotype']
      gt_pos = row['Genotype position']
      length = row['Length']
      cat = row['Category']
      if cat == 'ins':
        gapped_aligns.append(gt)
        continue
      if gt_pos == 'e':
        gapped_aligns.append('multiple deletion genotypes')
        continue

      gt_pos = int(gt_pos)
      gap_gt = gt[:cutsite - length + gt_pos] + '-'*length + gt[cutsite - length + gt_pos:]
      gapped_aligns.append(gap_gt)
    return gapped_aligns

  gap_gts = get_gapped_alignments(top10, stats)
  
  elements = []
  for gt, fq, length, cat in zip(gap_gts, fqs, lens, cats):
    elements.append(
      html.P(
        '%s %.1f %s %s' % (gt, fq, length, cat)
      )
    )

  return elements

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