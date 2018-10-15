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
        id = 'G_hidden-pred-df-stats',
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
      # Submit button
      ###################################################
      # Submit button
      html.Div([
        html.Button(
          'PREDICT REPAIR',
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
                    {'label': 'Cutsite', 'value': 'Cutsite'},
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
                  value = ['Precision', 'Frameshift (%)', 'MH strength', 'M.F. gt (%)'],
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
              )
            ], style = dict(
                transform = 'translateX(20px)',
                height = 60,
              )
            ),

          ],
        ),

        ##
        ], 
        style = dict(
          transform = 'translateX(90px)',
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
    width = '1150px',
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
# AWS S3 download callback
##
@app.callback(
  Output('G_hidden-pred-df-stats', 'children'),
  [Input('G_submit_button', 'n_clicks')]
)
def update_df_stats(nclicks):
  query_fn = 'hg38_0_mESC_SpCas9_0.csv'
  local_dir = 'local-s3/'
  s3.Bucket('indelphi-storage').download_file(query_fn, local_dir + query_fn)

  all_stats = pd.read_csv(local_dir + query_fn, index_col = 0)
  all_stats['ID'] = all_stats.index + 1
  all_stats['MH strength'] = np.log(all_stats['Phi'])
  return all_stats.iloc[:1000].to_csv()
  # return all_stats.iloc[:100].to_csv()
  # return all_stats.to_csv()

##
# Column selection and sorting callbacks
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
  Output('G_dropdown-columns', 'options'),
  [Input('G_hidden-pred-df-stats', 'children')],
  [State('G_dropdown-columns', 'options')]
  )
def update_columns_options(all_stats_string, prev_options):
  stats = pd.read_csv(StringIO(all_stats_string), index_col = 0)
  options = prev_options

  for d in ['Repairs to spec.', 'Deletes spec.', 'Dist. to POI']:
    td = {'label': d, 'value': d}
    if d in stats.columns:
      if td not in options:
        options.append(td)
    else:
      if td in options:
        options.remove(td)
  return options

@app.callback(
  Output('G_dropdown-columns', 'value'),
  [Input('G_dropdown-columns', 'options')],
  [State('G_dropdown-columns', 'value'),
   State('G_url', 'pathname'),
   State('G_row_dropdown-columns', 'n_clicks')]
  )
def update_columns_value(options, prev_value, url, n_clicks):
  value = prev_value
  all_options = [s['value'] for s in options]

  for td in ['Repairs to spec.', 'Deletes spec.', 'Dist. to POI']:
    if td in all_options:
      if td not in value:
        value.append(td)
    else:
      if td in value:
        value.remove(td)

  if n_clicks is None or n_clicks == 0:
    valid_flag, dd = lib.parse_valid_url_path_batch(url)
    if valid_flag:
      value = []
      alphabetical_options = sorted(all_options)
      for idx, flag in enumerate(dd['chosen_columns']):
        if flag == '1':
          value.append(alphabetical_options[idx])

  return value

##
# Stats table callbacks
## 
@app.callback(
  Output('G_table-stats', 'rows'), 
  [Input('G_hidden-pred-df-stats', 'children'),
   Input('G_dropdown-columns', 'value'),
   Input('G_dropdown-sortcol', 'value'),
   Input('G_sortdirection', 'value'),
  ])
def update_stats_table(all_stats_string, chosen_columns, sort_col, sort_direction):
  if all_stats_string == 'init':
    assert False, 'init'
  stats = pd.read_csv(StringIO(all_stats_string), index_col = 0)

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
    'Cas9 type',
    'Celltype',
    'Chromosome',
    'Cutsite distance to 3p boundary',
    'Cutsite distance to 5p boundary',
    'Exon end',
    'Exon number',
    'Exon start',
    'Exon strand',
    'Gene symbol',
    'Genome',
    'Local context',
    'Local cutsite',
    'gRNA',
    'gRNA strand w.r.t. exon strand',
    'kgID',
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
  return stats.to_dict('records')

@app.callback(
  Output('G_table-stats', 'selected_row_indices'),
  [Input('G_hidden-clickData', 'children'),
   Input('G_hidden-cache-submit-button', 'children'),
   Input('G_dropdown-columns', 'value'),
   Input('G_dropdown-sortcol', 'value'),
   Input('G_table-stats', 'rows')],
  [State('G_table-stats', 'selected_row_indices'),
   State('G_hidden-sort-module-interaction', 'children'),
   State('G_hidden-selected-id', 'children'),
   State('G_url', 'pathname'),
   State('G_postcomputation_settings', 'n_clicks'),
   State('G_plot-stats-div', 'n_clicks'),
   State('G_submit_button', 'n_clicks'),
   ])
def update_statstable_selected(clickData, submit_time, col_values, sortcol_value, rows, selected_row_indices, sort_time, prev_id, url, nc1, nc2, nc_submit):
  if not bool(nc1 or nc2) and nc_submit == 1:
    # On page load, select row from URL
    valid_flag, dd = lib.parse_valid_url_path_batch(url)
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
    df = pd.DataFrame(rows)
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
  [State('G_table-stats', 'rows')])
def update_hidden_selected_id(selected_idx, rows):
  if len(selected_idx) == 0:
    return ''
  idx = selected_idx[0]
  df = pd.DataFrame(rows)
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
    [Input('G_table-stats', 'rows'),
     Input('G_table-stats', 'selected_row_indices')])
def update_stats_plot(rows, selected_row_indices):
  df = pd.DataFrame(rows)
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
  fixedwidth_ids = lib.get_fixedwidth_ID(df['ID'])
  for idx, row in df.iterrows():
    row_text = '%s %s' % (row['gRNA'], fixedwidth_ids[idx])
    # row_text = '%s %s %s <a href="%s">details</a> %s' % (row['gRNA'], row['PAM'], row['gRNA orientation'], row['URL'], fixedwidth_ids[idx])
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
  fig['layout']['width'] = 275 + len(stats_cols) * 150
  fig['layout']['height'] = 150 + len(df) * 11
  fig['layout']['margin'] = {
    'l': 250,
    'r': 25,
    't': 0,
    'b': 150,
  }
  return fig

@app.callback(
    Output('G_hist-stats', 'figure'),
    [Input('G_table-stats', 'rows'),
     Input('G_table-stats', 'selected_row_indices')])
def update_hist_plot(rows, selected_row_indices):
  df = pd.DataFrame(rows)

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
  fig['layout']['width'] = 275 + len(stats_cols) * 150
  fig['layout']['height'] = 100
  fig['layout']['margin'] = {
    'l': 250,
    'r': 25,
    't': 0,
    # 't': 60,
    # 'b': 25,
    'b': 40,
  }
  return fig

