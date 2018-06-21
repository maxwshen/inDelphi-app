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
import plotly

import inDelphi
import generalStats
import lib

from indelphi_app import app

# init
inDelphi.init_model()
if not os.path.isdir('user-csvs/'):
  os.mkdir('user-csvs/')
else:
  subprocess.check_output('rm -rf user-csvs/*', shell = True)

# Remove these plotly modebar buttons to limit interactivity
modebarbuttons_2d = ['zoom2d', 'pan2d', 'select2d', 'lasso2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'resetScale2d', 'hoverClosestCartesian', 'hoverCompareCartesian', 'toggleSpikelines']

## Parameters
headerHeight = 200

###################################################################
###################################################################
##
# App layout
##
layout = html.Div([
  html.Div([

    ##
    # Hidden divs for light data storage
    ##
    html.Div(
      [
        html.Div(
          id = 'B2_hidden-pred-df-stats',
          children = 'init'
        ),
        html.Div(
          id = 'B2_hidden-cache-submit-button',
          children = '%s' % (time.time())
        ),

        # Datatable
        dt.DataTable(
          id = 'B2_table-stats',
          rows = [{}], # init rows
          selected_row_indices = [],
        ),

        dcc.Location(
          id = 'B2_url',
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
          'inDelphi batch mode dev',
          style = dict(
            textAlign = 'center',
          ),
        ),

        ###################################################
        # Sequence box
        ###################################################
        html.Div([
          dcc.Textarea(
            id = 'B2_textarea', 
            value = 'CTAGGGATGTGGCTGCATGCTACGTTGACACACCTACACTGCTCGAAGTAAATATACGAAGCGCGCGGCCTGGCCGGAGCCGTTCCGCATCGTCACGTGTTCGTTTACTGTTAATTGGTGGCACATAAGCAATATCGTAGTCCGTCAAATTCAGCCCTGTTATCCCCGGCGTTATGTGTCAAATGGCGTAGAACTGGATTGACTGTTTGACGGTACCTGCTGATCGGTACGGTGACCGAGAATCTGTCGGGCTATGTCACTAATACTTT',
            minLength = 70,  
            maxLength = 2000,  
            style = dict(
              fontFamily = 'monospace',
              fontSize = 16,
              resize = 'none',
              height = '60px',
              width = '800px',
            ),
          )],
          style = dict(
            verticalAlign = 'center',
            whiteSpace = 'nowrap',
            overflowX = 'auto',
            textAlign = 'center',
          ),
        ),

        ###################################################
        # PAM box
        ###################################################
        html.Div(
          [
            html.Strong(
              'Cas9 PAM: '
            ),
            dcc.Input(
              id = 'B2_textbox_pam', 
              size = 5,
              value = 'NGG',
              type = 'text',
              minlength = 2,
              maxlength = 6,
              autofocus = True,
              style = dict(
                fontFamily = 'monospace',
                fontSize = 14,
                height = '20px',
                width = '70px',
                verticalAlign = 'text-bottom',
              ),
            ),
            html.Strong(' ' * 3),
            html.Div(
              [
                html.Img(
                  src = '/staticfiles/tooltip_logo',
                  className = 'tooltiplogo',
                ),
                html.Span(
                  'Cutsite assumed 3nt upstream of PAM match. PAM must be 2-6 bp long. Supports IUPAC DNA encoding, ex: NNNRRT, NGG.',
                  className = 'tooltiptext'
                ),
              ], 
              className = 'tooltip',
            ),
          ],
          style = dict(
            textAlign = 'center',
          ),
        ),

        ###################################################
        # Click to run button
        ###################################################
        html.Div([
          html.Button(
            'Submit',
            id = 'B2_submit_button',
            style = dict(
            ),
          )],
          style = dict(
            textAlign = 'center',
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
        # Module: Detailed genotypes
        ###################################################

        # Hists
        html.Div(
          id = 'B2_hist-stats'
        ),

        # Plots
        html.Div(
          id = 'B2_plot-stats'
        ),

        # Sharable link, move this later
        html.Div(
          html.A(
            '🔗 Shareable link to page before computation', 
            id = 'B2_page-link'
          ),
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
# URL and header callbacks
##
@app.callback(
  Output('B2_textarea', 'value'),
  [Input('B2_url', 'pathname')],
  [State('B2_textarea', 'value')])
def update_textarea_from_url(url, default_value):
  valid_flag, textarea, pam = lib.parse_valid_url_path_batch(url)
  if valid_flag:
    return textarea
  return default_value

@app.callback(
  Output('B2_textbox_pam', 'value'),
  [Input('B2_url', 'pathname')],
  [State('B2_textbox_pam', 'value')])
def update_textarea_from_url(url, default_value):
  valid_flag, textarea, pam = lib.parse_valid_url_path_batch(url)
  if valid_flag:
    return pam
  return default_value

##
# Prediction callback
##
@app.callback(
  Output('B2_hidden-pred-df-stats', 'children'),
  [Input('B2_submit_button', 'n_clicks')],
  [State('B2_textarea', 'value'),
   State('B2_textbox_pam', 'value')])
def update_pred_df_stats(nclicks, seq, pam):
  dd = defaultdict(list)
  all_stats = pd.DataFrame()

  assert pam.count('N') != len(pam)
  assert 2 <= len(pam) <= 6

  # Search for gRNAs matching PAM
  seqs = [seq, lib.revcomp(seq)]
  cutsites = range(30, len(seq) - 30)
  for local_seq, grna_orient in zip(seqs, ['+', '-']):
    for local_cutsite in cutsites:
      cand_pam = local_seq[local_cutsite + 3 : local_cutsite + 3 + len(pam)]
      if lib.match(pam, cand_pam):
        dd['gRNA orientation'].append(grna_orient)
        dd['gRNA'].append(local_seq[local_cutsite - 17 : local_cutsite + 3])
        dd['PAM'].append(cand_pam)
        if grna_orient == '+':
          dd['Cutsite'].append(local_cutsite)
        else:
          dd['Cutsite'].append(len(seq) - local_cutsite)

        pred_df, stats = inDelphi.predict(local_seq, local_cutsite)
        all_stats = all_stats.append(stats, ignore_index = True)


  all_stats['Log phi'] = np.log(all_stats['Phi'])  
  # Drop
  drop_cols = [
    'Phi',
  ]
  all_stats = all_stats.drop(drop_cols, axis = 1)

  all_stats['ID'] = all_stats.index

  for col in dd:
    all_stats[col] = dd[col]
  return all_stats.to_csv()
  # return (pred_df.to_csv(), pd.DataFrame(stats, index = [0]).to_csv())


##
# Stats table callbacks
## 
@app.callback(
  Output('B2_table-stats', 'rows'), 
  [Input('B2_hidden-pred-df-stats', 'children')
  ])
def update_stats_table(all_stats_string):
  stats = pd.read_csv(StringIO(all_stats_string), index_col = 0)

  # Drop
  drop_cols = [
    'Reference sequence',
    '1-bp ins frequency',
    'MH del frequency',
    'MHless del frequency',
  ]
  stats = stats.drop(drop_cols, axis = 1)

  # Reformat floats
  stats_cols = list(stats.columns)
  nonstat_cols = ['ID', 'gRNA', 'gRNA orientation', 'PAM', 'Cutsite']
  for nonstat_col in nonstat_cols:
    stats_cols.remove(nonstat_col)
  for stat_col in stats_cols:
    if stat_col in ['Precision', 'Log phi']:
      stats[stat_col] = [float('%.2f' % (s)) for s in stats[stat_col]]    
    else:
      stats[stat_col] = [float('%.1f' % (s)) for s in stats[stat_col]]    

  # Reorder columns
  stats = stats[nonstat_cols + sorted(stats_cols)]

  return stats.to_dict('records')

@app.callback(
  Output('B2_table-stats', 'selected_row_indices'),
  [Input('B2_plot-stats-child', 'clickData')],
  [State('B2_table-stats', 'selected_row_indices')]
  )
def update_statstable_selected(clickData, selected_row_indices):
  # Update selections in table based on clicking plot
  if clickData:
    for point in clickData['points']:
      # Only allow selecting one point in plot-stats
      selected_row_indices = [point['pointNumber']]
  # Need to add: if hitting submit button, clear the selected rows. Otherwise, selecting a row M > number of rows N in new query, will fail
  return selected_row_indices

##
# Plot stats callback
##
@app.callback(
    Output('B2_plot-stats', 'children'),
    [Input('B2_table-stats', 'rows'),
     Input('B2_table-stats', 'selected_row_indices')])
def update_stats_plot(rows, selected_row_indices):
  try:
    df = pd.DataFrame(rows)
  except:
    # On page load, hide empty dash figure
    # Can put estimated loading time here
    return ''

  # Determine statistics to plot
  stats_cols = list(df.columns)
  nonstat_cols = ['ID', 'gRNA', 'gRNA orientation', 'PAM', 'Cutsite']
  for nonstat_col in nonstat_cols:
    stats_cols.remove(nonstat_col)
  stats_cols = sorted(stats_cols)

  fig = plotly.tools.make_subplots(
    rows = 1, cols = len(stats_cols),
    shared_yaxes = True)

  # Color selected markers

  if len(selected_row_indices) > 0:
    selected_row_index = selected_row_indices[0]
  else:
    selected_row_index = None
  selected_line = dict()

  # Generate each plot
  for idx, stats_col in enumerate(stats_cols):
    subplot_num = idx + 1
    marker = {'color': [lib.get_color(stats_col)] * len(df)}
    for i in (selected_row_indices or []):
      marker['color'][i] = '#000000'
    # Scatter
    fig.append_trace(
      go.Scattergl(
        x = df[stats_col],
        y = np.array(df.index) + 1,
        mode = 'markers',
        marker = marker,
        name = '',
      ), 
      1, subplot_num
    )
    if selected_row_index is not None:
      selected_line[subplot_num] = (df.index[selected_row_index], df[stats_col][selected_row_index])

  # Subplot formatting
  fig['layout']['yaxis1'].update(
    fixedrange = True,
    tickvals = np.arange(len(df.index)) + 1,
    ticktext = [str(s) for s in df['ID']],
    zeroline = True,
    zerolinewidth = 2,
    autorange = 'reversed',
    titlefont = dict(
      size = 10,
    ),
    range = [0, len(df)],
  )

  for idx, stats_col in enumerate(stats_cols):
    subplot_num = idx + 1
    [xmin, xmax] = lib.get_batch_statcol_xrange(df[stats_col], stats_col)
    fig['layout']['xaxis%s' % (subplot_num)].update(
      title = stats_col,
      fixedrange = True,
      showgrid = True,
      zeroline = False,
      titlefont = dict(
        size = 12,
      ),
      range = [xmin, xmax],
    )

    if selected_row_index is not None:
      fig['layout']['shapes'].append(
        lib.get_batch_select_line(
          x0 = selected_line[subplot_num][1],
          x1 = selected_line[subplot_num][1],
          y0 = 0,
          y1 = len(df),
          xref = 'x%s' % (subplot_num),
          yref = 'y1',
        )
      )
      fig['layout']['shapes'].append(
        lib.get_batch_select_line(
          x0 = xmin,
          x1 = xmax,
          y0 = selected_line[subplot_num][0] + 1,
          y1 = selected_line[subplot_num][0] + 1,
          xref = 'x%s' % (subplot_num),
          yref = 'y1',
        )
      )

  # Global figure formatting
  fig['layout']['showlegend'] = False
  fig['layout']['hovermode'] = 'y'
  fig['layout']['width'] = 150 * len(stats_cols)
  fig['layout']['height'] = 150 + len(df) * 11
  fig['layout']['margin'] = {
    'l': 25,
    'r': 25,
    't': 0,
    'b': 150,
  }
  child = dcc.Graph(
    id = 'B2_plot-stats-child',
    figure = fig,
    config = dict(
      modeBarButtonsToRemove = modebarbuttons_2d,
      displaylogo = False,
    ),
  )
  return child

@app.callback(
    Output('B2_hist-stats', 'children'),
    [Input('B2_table-stats', 'rows'),
     Input('B2_table-stats', 'selected_row_indices')])
def update_hist_plot(rows, selected_row_indices):
  try:
    df = pd.DataFrame(rows)
  except:
    # On page load, hide empty dash figure
    # Can put estimated loading time here
    return ''

  # Determine statistics to plot
  stats_cols = list(df.columns)
  nonstat_cols = ['ID', 'gRNA', 'gRNA orientation', 'PAM', 'Cutsite']
  for nonstat_col in nonstat_cols:
    stats_cols.remove(nonstat_col)
  stats_cols = sorted(stats_cols)

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
      selected_line[subplot_num] = (df.index[selected_row_index], df[stats_col][selected_row_index])

  # Subplot formatting

  for idx, stats_col in enumerate(stats_cols):
    subplot_num = idx + 1
    fig['layout']['yaxis%s' % (subplot_num)].update(
      fixedrange = True,
      showticklabels = False,
      showgrid = False,
    )
    fig['layout']['xaxis%s' % (subplot_num)].update(
      fixedrange = True,
      showgrid = True,
      zeroline = False,
      ticks = 'outside',
      ticklen = 3,
      range = lib.get_batch_statcol_xrange(df[stats_col], stats_col),
    )

    if selected_row_index is not None:
      fig['layout']['shapes'].append(
        lib.get_batch_select_line(
          x0 = selected_line[subplot_num][1],
          x1 = selected_line[subplot_num][1],
          y0 = 0,
          y1 = len(df) / 2.5,
          xref = 'x%s' % (subplot_num),
          yref = 'y1',
        )
      )

  # Global figure formatting
  fig['layout']['showlegend'] = False
  fig['layout']['width'] = 150 * len(stats_cols)
  fig['layout']['height'] = 150
  fig['layout']['margin'] = {
    'l': 25,
    'r': 25,
    't': 60,
    'b': 25,
  }
  child = dcc.Graph(
    id = 'B2_hist-stats-child',
    figure = fig,
    config = dict(
      modeBarButtonsToRemove = modebarbuttons_2d,
      displaylogo = False,
    ),
  )
  return child



##
# Page link callback
##
@app.callback(
  Output('B2_page-link', 'href'),
  [Input('B2_textarea', 'value'),
   Input('B2_textbox_pam', 'value'),
  ])
def update_pagelink(textarea, pam):
  return 'https://dev.crisprindelphi.design%s' % (lib.encode_dna_to_url_path_batch(textarea, pam))