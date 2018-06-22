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

        # Datatable
        dt.DataTable(
          id = 'B_table-stats',
          rows = [{}], # init rows
          selected_row_indices = [],
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
            id = 'B_textarea', 
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
              id = 'B_textbox_pam', 
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

        # Multi drop down to select columns
        dcc.Dropdown(
          id = 'B_dropdown-columns',
          options = [
            {'label': 'Cutsite', 'value': 'Cutsite'},
            {'label': 'Precision', 'value': 'Precision'},
            {'label': 'Frameshift (%)', 'value': 'Frameshift (%)'},
            {'label': 'Frame +0 (%)', 'value': 'Frame +0 (%)'},
            {'label': 'Frame +1 (%)', 'value': 'Frame +1 (%)'},
            {'label': 'Frame +2 (%)', 'value': 'Frame +2 (%)'},
            {'label': 'Log phi (microhomology strength)', 'value': 'Log phi'},
            {'label': 'Most frequent genotype (%)', 'value': 'M.F. gt (%)'},
            {'label': 'Most frequent deletion (%)', 'value': 'M.F. del (%)'},
            {'label': 'Most frequent insertion (%)', 'value': 'M.F. ins (%)'},
            {'label': 'Expected indel length', 'value': 'Exp. indel len'},
          ],
          multi = True,
          searchable = False,
          clearable = False,
          value = ['Cutsite', 'Precision', 'Frameshift (%)', 'Log phi', 'M.F. gt (%)']
        ),

        # Sorting columns
        dcc.Dropdown(
          id = 'B_dropdown-sortcol',
          options = [],
          searchable = False,
          clearable = False,
        ),
        # Sort direction
        dcc.RadioItems(
          id = 'B_sortdirection',
          options = [
            {'label': 'Up', 'value': 'Up'},
            {'label': 'Down', 'value': 'Down'},
          ],
          labelStyle = {'display': 'inline-block'}
        ),

        # Hists
        html.Div(
          id = 'B_hist-stats'
        ),

        # Plots
        html.Div(
          id = 'B_plot-stats'
        ),

        # Sharable link, move this later
        html.Div(
          html.A(
            '🔗 Shareable link to page before computation', 
            id = 'B_page-link'
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
  Output('B_textarea', 'value'),
  [Input('B_url', 'pathname')],
  [State('B_textarea', 'value')])
def update_textarea_from_url(url, default_value):
  valid_flag, textarea, pam = lib.parse_valid_url_path_batch(url)
  if valid_flag:
    return textarea
  return default_value

@app.callback(
  Output('B_textbox_pam', 'value'),
  [Input('B_url', 'pathname')],
  [State('B_textbox_pam', 'value')])
def update_pam_from_url(url, default_value):
  valid_flag, textarea, pam = lib.parse_valid_url_path_batch(url)
  if valid_flag:
    return pam
  return default_value

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

  all_stats['ID'] = all_stats.index + 1

  for col in dd:
    all_stats[col] = dd[col]
  return all_stats.to_csv()
  # return (pred_df.to_csv(), pd.DataFrame(stats, index = [0]).to_csv())

##
# Column selection and sorting callbacks
##
@app.callback(
  Output('B_dropdown-sortcol', 'options'),
  [Input('B_dropdown-columns', 'value')])
def update_sortcol_options(values):
  options = []
  for value in values:
    options.append({'label': value, 'value': value})
  return options

##
# Stats table callbacks
## 
@app.callback(
  Output('B_table-stats', 'rows'), 
  [Input('B_hidden-pred-df-stats', 'children'),
   Input('B_dropdown-columns', 'value'),
   Input('B_dropdown-sortcol', 'value'),
   Input('B_sortdirection', 'value'),
  ])
def update_stats_table(all_stats_string, chosen_columns, sort_col, sort_direction):
  stats = pd.read_csv(StringIO(all_stats_string), index_col = 0)

  # Drop extra cols
  drop_cols = [
    'Reference sequence',
    '1-bp ins frequency',
    'MH del frequency',
    'MHless del frequency',
  ]
  stats = stats.drop(drop_cols, axis = 1)

  # Rename to shorter versions
  stats = lib.rename_batch_columns(stats)

  # Sort by, if possible
  if sort_col is not None and sort_direction is not None:
    if sort_direction == 'Up':
      ascending_flag = True
    else:
      ascending_flag = False
    print(stats.columns)
    stats = stats.sort_values(by = sort_col, ascending = ascending_flag)

  # Reformat floats
  stats_cols = list(stats.columns)
  nonstat_cols = ['ID', 'gRNA', 'gRNA orientation', 'PAM']
  for nonstat_col in nonstat_cols:
    stats_cols.remove(nonstat_col)
  for stat_col in stats_cols:
    # Filter down to selected columns
    if stat_col not in chosen_columns:
      stats.drop(stat_col, axis = 1, inplace = True)
      continue
    # Reformat
    if stat_col in ['Precision', 'Log phi']:
      stats[stat_col] = [float('%.2f' % (s)) for s in stats[stat_col]]    
    else:
      stats[stat_col] = [float('%.1f' % (s)) for s in stats[stat_col]]    

  # Reorder columns
  stats = stats[nonstat_cols + lib.order_chosen_columns(chosen_columns)]

  return stats.to_dict('records')

@app.callback(
  Output('B_table-stats', 'selected_row_indices'),
  [Input('B_plot-stats-child', 'clickData')],
  [State('B_table-stats', 'selected_row_indices')]
  )
def update_statstable_selected(clickData, selected_row_indices):
  # Update selections in table based on clicking plot
  if clickData:
    clicked_idx = clickData['points'][0]['pointNumber']
    if selected_row_indices != [clicked_idx]:
      selected_row_indices = [clicked_idx]
    else:
      # Point already selected, user clicked on same point twice:
      # so, deselect
      selected_row_indices = []
      # Only allow selecting one point in plot-stats
  # Need to add: if hitting submit button, clear the selected rows. Otherwise, selecting a row M > number of rows N in new query, will fail
  # Need to add: If changing sort col or direction, clear the selected rows. Otherwise, the wrong row is selected after sorting. Preferably, keep the selected row and update the index.
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
    # Can put estimated loading time here
    return ''

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
    tickvals = np.arange(1, len(df.index) + 1),
    ticktext = [str(s) for s in df['ID']],
    zeroline = True,
    zerolinewidth = 2,
    autorange = 'reversed',
    titlefont = dict(
      size = 10,
    ),
    range = [0, len(df)],
  )

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
  # fig['layout']['spikedistance'] = -1
  fig['layout']['width'] = 50 + len(stats_cols) * 150
  fig['layout']['height'] = 150 + len(df) * 11
  fig['layout']['margin'] = {
    'l': 25,
    'r': 25,
    't': 0,
    'b': 150,
  }
  child = dcc.Graph(
    id = 'B_plot-stats-child',
    figure = fig,
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
def update_hist_plot(rows, selected_row_indices):
  try:
    df = pd.DataFrame(rows)
  except:
    # On page load, hide empty dash figure
    # Can put estimated loading time here
    return ''

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
      selected_line[subplot_num] = (df.index[selected_row_index], df[stats_col][selected_row_index])

  # Subplot formatting

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
  fig['layout']['width'] = 50 + len(stats_cols) * 150
  fig['layout']['height'] = 140
  fig['layout']['margin'] = {
    'l': 25,
    'r': 25,
    't': 60,
    # 'b': 25,
    'b': 40,
  }
  child = dcc.Graph(
    id = 'B_hist-stats-child',
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
  Output('B_page-link', 'href'),
  [Input('B_textarea', 'value'),
   Input('B_textbox_pam', 'value'),
  ])
def update_pagelink(textarea, pam):
  return 'https://dev.crisprindelphi.design%s' % (lib.encode_dna_to_url_path_batch(textarea, pam))