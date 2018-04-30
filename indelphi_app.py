import pickle, copy, os, datetime, subprocess
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy.stats import entropy

import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_table_experiments as dt
from dash.dependencies import Input, Output, State
import plotly.graph_objs as go
import flask

import inDelphi
import generalStats


# Dash server setup
app = dash.Dash('')
server = app.server

# init
inDelphi.init_model()
generalStats.init_all_traces()
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
                size = 35,
                type = 'text',
                value = 'TAGTTTCTAGCACAGCGCTGGTGTGGC', 
                autofocus = True,
                style = dict(
                  textAlign = 'right',
                  fontFamily = 'monospace',
                  fontSize = 16,
                ),
              )
            ],
            style = dict(
              marginLeft = '55px',
            ),
            className = 'six columns',
          ),

          # Right box
          html.Div(
            [
              dcc.Input(
                id = 'textbox2', 
                size = 35,
                type = 'text',
                value = 'GTGTGGCTGAAGGCATAGTAATTCTGA', 
                style = dict(
                  textAlign = 'left',
                  fontFamily = 'monospace',
                  fontSize = 16,
                ),
              ),
            ],
            style = dict(
              marginLeft = '-40px',
            ),
            className = 'six columns',
          ),
        ], 
        className = 'row'),

        html.Div(
          'DSB',
          style = dict(
            textAlign = 'center',
          )
        ),
      ],
      style = dict(
        position = 'fixed',
        backgroundColor = 'white',
        borderBottom = '5px solid #DDDDDD',
        zIndex = 1e6,
        width = '800px',
        margin = '0 auto',
        height = headerHeight,
        marginTop = '-%spx' % (headerHeight),
      ),
    ),

    ##
    # Body / plots
    ##
    html.Div(
      [
        # First row
        html.Div(
          [
            # Frameshift
            html.Div(
              [
                dcc.Graph(
                  id = 'plot-fs',
                  style = dict(
                    height = 300, 
                    width = 250,
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
                    width = 600,
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

        dcc.Graph(
          id = 'plot-genstats-precision',
          style = dict(
            height = 300, 
            width = 500,
          ),
          config = dict(
            modeBarButtonsToRemove = modebarbuttons_2d,
            displaylogo = False,
          ),
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

        html.A(
          'Download CSV', 
          id = 'csv-download-link'
        ),

        html.Div(
          'Copyright MIT 2018.\nAll Rights Reserved.',
          style = dict(
            textAlign = 'center',
          )
        ),
      ],
      # body style
      style = dict(
        marginTop = '%spx' % (headerHeight),
      ),
    ),
    ##
  ], 
    style = dict(
      width = '800px',
      margin = '0 auto',
    )
  ),
], # Shadow
style = dict(
  width = '900px',
  margin = '0 auto',
  boxShadow = '0 4px 8px 0 rgba(0, 0, 0, 0.2), 0 6px 20px 0 rgba(0, 0, 0, 0.19)',
)
)

#######################################################################
#########################      CALLBACKS      #########################
#######################################################################
##
# Header Callbacks
##


##
# General stats callbacks
##
@app.callback(
  Output('plot-genstats-precision', 'figure'),
  [Input('textbox1', 'value'),
   Input('textbox2', 'value'),
  ])
def cb_plot_genstats_precision(text1, text2):
  seq = text1 + text2
  cutsite = len(text1)
  pred_df, stats = inDelphi.predict(seq, cutsite)
  xval = inDelphi.get_precision(pred_df)
  return dict(
    data = [
      generalStats.trace_precision,
    ],
    layout = generalStats.layout('Precision score', xval),
  )

##
# Indel length and frameshift callbacks
@app.callback(
  Output('plot-indel-len', 'figure'),
  [Input('textbox1', 'value'),
   Input('textbox2', 'value'),
  ])
def cb_plot_indel_len(text1, text2):
  seq = text1 + text2
  cutsite = len(text1)
  pred_df, stats = inDelphi.predict(seq, cutsite)
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
        t = 10,
      ),
    ),
  )

@app.callback(
  Output('plot-fs', 'figure'),
  [Input('textbox1', 'value'),
   Input('textbox2', 'value'),
  ])
def cb_plot_fs(text1, text2):
  seq = text1 + text2
  cutsite = len(text1)
  pred_df, stats = inDelphi.predict(seq, cutsite)
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
        t = 10,
      ),
    ),
  )



##
# Genotype table callbacks
## 
@app.callback(
  Output('table-genotypes', 'rows'), 
  [Input('textbox1', 'value'),
   Input('textbox2', 'value'),
  ])
def cb_update_genotype_table(text1, text2):
  seq = text1 + text2
  cutsite = len(text1)
  pred_df, stats = inDelphi.predict(seq, cutsite)
  inDelphi.add_genotype_column(pred_df, stats)
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
        x = df['Genotype'],
        y = df['Predicted frequency'],
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
  [Input('textbox1', 'value'),
   Input('textbox2', 'value'),
  ])
def update_link(text1, text2):
  seq = text1 + text2
  cutsite = len(text1)

  pred_df, stats = inDelphi.predict(seq, cutsite)
  inDelphi.add_genotype_column(pred_df, stats)

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

###################################################################
###################################################################
# CSS
app.css.append_css({'external_url': 'https://codepen.io/chriddyp/pen/bWLwgP.css'})

if __name__ == '__main__':
  app.run_server()