import dash
import dash_core_components as dcc
import dash_html_components as html

# Dash server setup
app = dash.Dash('')
server = app.server

##
# App layout
##
app.layout = html.Div([
  html.H1('inDelphi'),
  html.Div('Description'),

  # Input text box 1
  html.Label('DNA sequence 1'),
  dcc.Input(id = 'textbox1', value = 'ACGT', type = 'text'),

  # Input text box 2
  html.Label('DNA sequence 2'),
  dcc.Input(id = 'textbox2', value = 'ACGT', type = 'text'),

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


# CSS
app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"})

if __name__ == '__main__':
  app.run_server()