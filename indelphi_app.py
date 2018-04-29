import pickle, copy, os, datetime
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

# Remove these plotly modebar buttons to limit interactivity
modebarbuttons_2d = ['zoom2d', 'pan2d', 'select2d', 'lasso2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'resetScale2d', 'hoverClosestCartesian', 'hoverCompareCartesian', 'toggleSpikelines']



###################################################################
###################################################################
##
# App layout
##
app.layout = html.Div([

  # Non-scrolling header
  html.Div(
    [
      html.H1('inDelphi'),
      html.Div('Description',
        style = dict(
          textAlign = 'center',
        )
      ),

      # Input text box 1
      html.Div([
        html.Label('DNA sequence 1'),
        dcc.Input(id = 'textbox1', 
                  value = 'TAGTTTCTAGCACAGCGCTGGTGTGGC', 
                  type = 'text'),
        ], 
        style = dict(
          textAlign = 'center',
        )
      ),

      # Input text box 2
      html.Label('DNA sequence 2'),
      dcc.Input(id = 'textbox2', 
                value = 'GTGTGGCTGAAGGCATAGTAATTCTGA', 
                type = 'text'),

      html.Div(id = 'seq_display'),
    ],
    style = dict(
      # position = 'fixed',
    ),
  ),

  html.Table(id = 'fs_table'),

  dcc.Graph(
    id = 'plot-fs',
    style = dict(
      height = 400, 
      width = 350,
    ),
    config = dict(
      modeBarButtonsToRemove = modebarbuttons_2d,
      displaylogo = False,
    ),
  ),

  dcc.Graph(
    id = 'plot-genstats-precision',
    style = dict(
      height = 400, 
      width = 500,
    ),
    config = dict(
      modeBarButtonsToRemove = modebarbuttons_2d,
      displaylogo = False,
    ),
  ),

  dcc.Graph(
    id = 'plot-indel-len',
    style = dict(
      height = 400, 
      width = 600,
    ),
    config = dict(
      modeBarButtonsToRemove = modebarbuttons_2d,
      displaylogo = False,
    ),
  ),

  html.Div('-----' * 5),

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

], 
  style = dict(
    width = '800px',
    margin = '0 auto',
  )
)

###################################################################
###################################################################
##
# Callbacks
##
@app.callback(
  Output('seq_display', 'children'),
  [Input('textbox1', 'value'),
   Input('textbox2', 'value'),
  ])
def cb_displaytext_seq(text1, text2):
  seq = text1 + text2
  return seq

@app.callback(
  Output('fs_table', 'children'),
  [Input('textbox1', 'value'),
   Input('textbox2', 'value'),
  ])
def cb_table_fs(text1, text2):
  seq = text1 + text2
  cutsite = len(text1)
  ans = inDelphi.predict(seq, cutsite)
  pred_df, stats = ans
  fs_df = inDelphi.get_frameshift_fqs(pred_df)
  max_rows = 10
  return html.Table(
    # Header
    [html.Tr([html.Th(col) for col in fs_df.columns])] +

    # Body
    [html.Tr([
        html.Td(fs_df.iloc[i][col]) for col in fs_df.columns
    ]) for i in range(min(len(fs_df), max_rows))]
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
        showgrid = False,
        zeroline = False,
        showline = True,
        ticks = 'outside',
        tick0 = 0,
        ticklen = 3,
        tickwidth = 0.5,
      ),
      font = dict(
        family = 'Arial',
      ),
    ),
  )


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

  X = [-1*int(s) for s in lendf['Indel length']]
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
    # layout = 
  )

@app.callback(
  Output('table-genotypes', 'rows'), 
  [Input('textbox1', 'value'),
   Input('textbox2', 'value'),
  ])
def update_datatable(text1, text2):
  seq = text1 + text2
  cutsite = len(text1)
  pred_df, stats = inDelphi.predict(seq, cutsite)
  return pred_df.to_dict('records')

@app.callback(
  Output('csv-download-link', 'href'), 
  [Input('textbox1', 'value'),
   Input('textbox2', 'value'),
  ])
def update_link(text1, text2):
  seq = text1 + text2
  cutsite = len(text1)
  pred_df, stats = inDelphi.predict(seq, cutsite)
  
  time = str(datetime.datetime.now()).replace(' ', '_').replace(':', '-')
  link_fn = '/dash/urlToDownload?value={}'.format(time)
  pred_df.to_csv('%s.csv' % (time))
  return link_fn

@app.server.route('/dash/urlToDownload') 
def download_csv():
  value = flask.request.args.get('value')
  # create a dynamic csv or file here using `StringIO` 
  # (instead of writing to the file system)
  local_csv_fn = value.split('/')[-1]
  return flask.send_file(
    # value,
    open(local_csv_fn + '.csv', 'rb'),
    mimetype = 'text/csv',
    attachment_filename = 'inDelphi_output.csv',
    as_attachment = True,
  )

###################################################################
###################################################################
# CSS
app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"})

if __name__ == '__main__':
  app.run_server()