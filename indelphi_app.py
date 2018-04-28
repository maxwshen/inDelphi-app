import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go

import inDelphi


# Dash server setup
app = dash.Dash('')
server = app.server

# Model init
inDelphi.init_model()

###################################################################
###################################################################
##
# App layout
##
app.layout = html.Div([
  html.H1('inDelphi'),
  html.Div('Description'),

  # Input text box 1
  html.Label('DNA sequence 1'),
  dcc.Input(id = 'textbox1', 
            value = 'TAGTTTCTAGCACAGCGCTGGTGTGGC', 
            type = 'text'),

  # Input text box 2
  html.Label('DNA sequence 2'),
  dcc.Input(id = 'textbox2', 
            value = 'GTGTGGCTGAAGGCATAGTAATTCTGA', 
            type = 'text'),

  html.Div(id = 'seq_display'),

  html.Table(id = 'fs_table'),

  dcc.Graph(
    style = {'height': 300},
    id = 'plot-fs'
  ),

  html.Div('-----' * 5),

])

###################################################################
###################################################################
##
# Callbacks
##
@app.callback(
  dash.dependencies.Output('seq_display', 'children'),
  [dash.dependencies.Input('textbox1', 'value'),
   dash.dependencies.Input('textbox2', 'value'),
  ])
def cb_displayseq(text1, text2):
  seq = text1 + text2
  return seq

@app.callback(
  dash.dependencies.Output('fs_table', 'children'),
  [dash.dependencies.Input('textbox1', 'value'),
   dash.dependencies.Input('textbox2', 'value'),
  ])
def cb_display_fstable(text1, text2):
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
  dash.dependencies.Output('plot-fs', 'figure'),
  [dash.dependencies.Input('textbox1', 'value'),
   dash.dependencies.Input('textbox2', 'value'),
  ])
def cb_display_plot_fs(text1, text2):
  seq = text1 + text2
  cutsite = len(text1)
  ans = inDelphi.predict(seq, cutsite)
  pred_df, stats = ans
  fs_df = inDelphi.get_frameshift_fqs(pred_df)
  X = ['+0', '+1', '+2']
  Y = [float(fs_df[fs_df['Frame'] == s]['Predicted frequency']) for s in X]
  return {
    'data': [
      go.Bar(
        x = X,
        y = Y,
        marker = go.Marker(color = 'rgb(55, 83, 109)'),
      ),
    ],
    'layout': go.Layout(
      title = 'Frameshift frequency',
      margin = go.Margin(l = 40, r = 0, t = 40, b = 30)
    ),
  }

###################################################################
###################################################################
# CSS
app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"})

if __name__ == '__main__':
  app.run_server()