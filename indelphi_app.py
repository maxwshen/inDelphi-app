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

# Remove these plotly modebar buttons to limit interactivity
modebarbuttons_2d = ['zoom2d', 'pan2d', 'select2d', 'lasso2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'resetScale2d', 'hoverClosestCartesian', 'hoverCompareCartesian', 'toggleSpikelines']

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
def cb_displaytext_seq(text1, text2):
  seq = text1 + text2
  return seq

@app.callback(
  dash.dependencies.Output('fs_table', 'children'),
  [dash.dependencies.Input('textbox1', 'value'),
   dash.dependencies.Input('textbox2', 'value'),
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
  dash.dependencies.Output('plot-fs', 'figure'),
  [dash.dependencies.Input('textbox1', 'value'),
   dash.dependencies.Input('textbox2', 'value'),
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
        textposition = 'auto',
        opacity = 0.6,
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

###################################################################
###################################################################
# CSS
app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"})

if __name__ == '__main__':
  app.run_server()