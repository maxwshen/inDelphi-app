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
          id = 'B_hidden-pred-df-stats',
          children = 'init'
        ),
        html.Div(
          id = 'B_hidden-cache-submit-button',
          children = '%s' % (time.time())
        ),

        dcc.Location(
          id = 'B_url',
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
          'inDelphi batch mode',
          style = dict(
            textAlign = 'center',
          ),
        ),

        ###################################################
        # Sequence box
        ###################################################
        html.Div([
          dcc.Textarea(
            id = 'B_textarea', 
            value = 'CTAGGGATGTGGCTGCATGCTACGTTGACACACCTACACTGCTCGAAGTAAATATACGAAGCGCGCGGCCTGGCCGGAGCCGTTCCGCATCGTCACGTGTTCGTTTACTGTTAATTGGTGGCACATAAGCAATATCGTAGTCCGTCAAATTCAGCCCTGTTATCCCCGGCGTTATGTGTCAAATGGCGTAGAACTGGATTGACTGTTTGACGGTACCTGCTGATCGGTACGGTGACCGAGAATCTGTCGGGCTATGTCACTAATACTTTCCAAACGCCCCGTATCGATGCTGAACGAATCGATGCACGCTCCCGTCTTTGAAAACGCATAAACATACAAGTGGACAGATGATGGGTACGGGCCTCTAATACATCCAACACTCTACGCCCTCTTCAAGAGCTAGAAGGGCACCCTGCAGTTGGAAAGGGAATTATTTCGTAAGGCGAGCCCATACCGTCATTCATGCGGAAGAGTTAACACGATTGGAAGTAGGAATAGTT',
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
        html.Div([
          dcc.Input(
            id = 'B_textbox_pam', 
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
              verticalAlign = 'top',
            ),
          )],
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
            id = 'B_submit_button',
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
  
        dt.DataTable(
          id = 'B_table-stats',
          rows = [{}], # init rows
          row_selectable = True,
          # row_single_select = True,
          filterable = True,
          sortable = True,
          selected_row_indices = [],
        ),

        html.Div(
          [
            # Plots: Scatter
            html.Div(
              id = 'B_plot-stats',
              className = 'nine columns',
            ),

            # Plots: Histograms
            html.Div(
              id = 'B_hist-stats',
              className = 'three columns',
            ),

          ],
          className = 'row',
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
# Prediction callback
##
@app.callback(
  Output('B_hidden-pred-df-stats', 'children'),
  [Input('B_submit_button', 'n_clicks')],
  [State('B_textarea', 'value'),
   State('B_textbox_pam', 'value')])
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
  Output('B_table-stats', 'rows'), 
  [Input('B_hidden-pred-df-stats', 'children')
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
  Output('B_table-stats', 'selected_row_indices'),
  [Input('B_plot-stats', 'clickData')],
  [State('B_table-stats', 'selected_row_indices')]
  )
def update_statstable_selected(clickData, selected_row_indices):
  # Update selections in table based on clicking plot
  if clickData:
    for point in clickData['points']:
      if point['pointNumber'] in selected_row_indices:
        selected_row_indices.remove(point['pointNumber'])
      else:
        selected_row_indices.append(point['pointNumber'])
  return selected_row_indices

##
# Plot stats callback
##
@app.callback(
    Output('B_plot-stats', 'children'),
    [Input('B_table-stats', 'rows'),
     Input('B_table-stats', 'selected_row_indices')])
def update_stats_plot(rows, selected_row_indices):
  try:
    df = pd.DataFrame(rows)
  except:
    # On page load, hide empty dash figure
    return ''

  # Determine statistics to plot
  stats_cols = list(df.columns)
  nonstat_cols = ['ID', 'gRNA', 'gRNA orientation', 'PAM', 'Cutsite']
  for nonstat_col in nonstat_cols:
    stats_cols.remove(nonstat_col)
  stats_cols = sorted(stats_cols)

  fig = plotly.tools.make_subplots(
    rows = len(stats_cols), cols = 1)

  # Color selected markers
  marker = {'color': ['#0074D9'] * len(df)}
  for i in (selected_row_indices or []):
    marker['color'][i] = '#FF851B'

  if len(selected_row_indices) > 0:
    selected_row_index = selected_row_indices[0]
  else:
    selected_row_index = None
  selected_line = dict()

  # Generate each plot
  for idx, stats_col in enumerate(stats_cols):
    subplot_num = idx + 1
    # Scatter
    fig.append_trace(
      go.Scatter(
        x = np.array(df.index) + 1,
        y = df[stats_col],
        mode = 'markers',
        marker = marker,
        name = '',
        xaxis = 'x_%s' % (subplot_num),
        yaxis = 'y_%s' % (subplot_num),
      ), 
      subplot_num, 1
    )
    if selected_row_index is not None:
      selected_line[subplot_num] = (df.index[selected_row_index], df[stats_col][selected_row_index])

  # Subplot formatting
  for idx, stats_col in enumerate(stats_cols):
    subplot_num = idx + 1
    fig['layout']['xaxis%s' % (subplot_num)].update(
      fixedrange = True,
      tickvals = np.arange(len(df.index)) + 1,
      ticktext = [str(s) for s in df['ID']],
      zeroline = True,
    )

    fig['layout']['yaxis%s' % (subplot_num)].update(
      title = stats_col,
      fixedrange = True,
      showgrid = True,
      zeroline = False,
    )

    if selected_row_index is not None:
      fig['layout']['shapes'].append(
        dict(
          type = 'line',
          xref = 'x%s' % (subplot_num),
          yref = 'y%s' % (subplot_num),
          x0 = 0,
          x1 = len(df),
          y0 = selected_line[subplot_num][1],
          y1 = selected_line[subplot_num][1],
          opacity = 0.8,
          line = dict(
            color = 'rgb(33, 33, 33)',
            width = 1,
          ),
        )
      )

  # Global figure formatting
  fig['layout']['showlegend'] = False
  fig['layout']['height'] = 200 * len(stats_cols)
  fig['layout']['margin'] = {
    'l': 40,
    'r': 10,
    't': 60,
    'b': 200
  }
  child = dcc.Graph(
    figure = fig,
    id = 'B_plot-stats-child',
    config = dict(
      modeBarButtonsToRemove = modebarbuttons_2d,
      displaylogo = False,
    ),
  )
  return child

@app.callback(
    Output('B_hist-stats', 'children'),
    [Input('B_table-stats', 'rows'),
     Input('B_table-stats', 'selected_row_indices')])
def update_stats_plot(rows, selected_row_indices):
  try:
    df = pd.DataFrame(rows)
  except:
    return ''

  # Determine statistics to plot
  stats_cols = list(df.columns)
  nonstat_cols = ['ID', 'gRNA', 'gRNA orientation', 'PAM', 'Cutsite']
  for nonstat_col in nonstat_cols:
    stats_cols.remove(nonstat_col)
  stats_cols = sorted(stats_cols)

  fig = plotly.tools.make_subplots(
    rows = len(stats_cols), cols = 1)

  if len(selected_row_indices) > 0:
    selected_row_index = selected_row_indices[0]
  else:
    selected_row_index = None
  selected_line = dict()

  # Generate each plot
  for idx, stats_col in enumerate(stats_cols):
    subplot_num = idx + 1
    # Hist
    fig.append_trace(
      go.Histogram(
        y = df[stats_col],
        # 'marker': marker,
        name = '',
        xaxis = 'x_%s' % (subplot_num),
        yaxis = 'y_%s' % (subplot_num),
      ), 
      subplot_num, 1
    )
    if selected_row_index is not None:
      selected_line[subplot_num] = (df.index[selected_row_index], df[stats_col][selected_row_index])

  # Subplot formatting
  for idx, stats_col in enumerate(stats_cols):
    subplot_num = idx + 1
    fig['layout']['xaxis%s' % (subplot_num)].update(
      fixedrange = True,
    )
    fig['layout']['yaxis%s' % (subplot_num)].update(
      title = stats_col,
      fixedrange = True,
      showgrid = True,
      zeroline = False,
    )

    # Draw horizontal line for selection
    if selected_row_index is not None:
      fig['layout']['shapes'].append(
        dict(
          type = 'line',
          xref = 'x%s' % (subplot_num),
          yref = 'y%s' % (subplot_num),
          x0 = 0,
          x1 = len(df) / 2,
          y0 = selected_line[subplot_num][1],
          y1 = selected_line[subplot_num][1],
          opacity = 0.8,
          line = dict(
            color = 'rgb(33, 33, 33)',
            width = 1,
          ),
        )
      )

  # Global figure formatting
  fig['layout']['showlegend'] = False
  fig['layout']['height'] = 200 * len(stats_cols)
  fig['layout']['margin'] = {
    'l': 40,
    'r': 10,
    't': 60,
    'b': 200
  }
  child = dcc.Graph(
    figure = fig,
    id = 'B_hist-stats-child',
    config = dict(
      modeBarButtonsToRemove = modebarbuttons_2d,
      displaylogo = False,
    ),
  )
  return child
