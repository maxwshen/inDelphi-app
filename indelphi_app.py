import dash
import dash_core_components as dcc
import dash_html_components as html
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

  html.Div('-----' * 5),

  # remainder
  dcc.RadioItems(
    id='dropdown-a',
    options=[{'label': i, 'value': i} for i in ['Canada', 'USA', 'Mexico']],
    value='Canada'
  ),
  html.Div(id='output-a'),

  dcc.RadioItems(
    id='dropdown-b',
    options=[{'label': i, 'value': i} for i in ['MTL', 'NYC', 'SF']],
    value='MTL'
  ),
  html.Div(id='output-b')

])

###################################################################
###################################################################
##
# Callbacks
##
@app.callback(
  dash.dependencies.Output('output-a', 'children'),
  [dash.dependencies.Input('dropdown-a', 'value')])
def callback_a(dropdown_value):
  return 'You\'ve selected "{}"'.format(dropdown_value)


@app.callback(
  dash.dependencies.Output('output-b', 'children'),
  [dash.dependencies.Input('dropdown-a', 'value')])
def callback_b(dropdown_value):
  return 'You\'ve selected "{}"'.format(dropdown_value)

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


###################################################################
###################################################################
# CSS
app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"})

if __name__ == '__main__':
  app.run_server()