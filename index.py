import os

import dash_core_components as dcc
import dash_html_components as html
import dash_table_experiments as dt
from dash.dependencies import Input, Output, State
import flask

from indelphi_app import app
from apps import app_single, app_batch

###################################################################
###################################################################
# Layout
app.layout = html.Div([
    html.Div(id = 'master-page-content'),
    dcc.Location(id = 'master-url', refresh = False),

    # Hidden datatable for proper rendering
    # https://community.plot.ly/t/display-tables-in-dash/4707/40?u=chriddyp
    html.Div(dt.DataTable(rows=[{}]), style={'display': 'none'})
])

app.title = 'inDelphi'


###################################################################
###################################################################
# Serve pages
@app.callback(
  Output('master-page-content', 'children'),
  [Input('master-url', 'pathname')]
)
def display_page(pathname):
  # return app_single.layout
  print(pathname)
  if pathname is None or pathname == '/':
    return app_single.layout
  elif pathname[:len('/single')] == '/single':
    return app_single.layout
  elif pathname[:len('/batch')] == '/batch':
    return app_batch.layout
  else:
    return app_single.layout
  #   # return '404'

###################################################################
###################################################################
# CSS
css_directory = os.getcwd()
@app.server.route('/static/<stylesheet>')
def serve_stylesheet(stylesheet):
  if stylesheet not in stylesheets:
    raise Exception(
      '"{}" is excluded from the allowed static files'.format(
        stylesheet
      )
    )
  return flask.send_from_directory(css_directory, stylesheet)

stylesheets = ['stylesheet.css']
for stylesheet in stylesheets:
  app.css.append_css({'external_url': '/static/{}'.format(stylesheet)})

# As of 0.22.0, served automatically from /assets/

# Favicon
@app.server.route('/favicon.ico')
def favicon():
  return flask.send_from_directory(os.get_cwd() + '/staticfiles/', 'favicon_CEL_icon.ico', mimetype = 'image/vnd.microsoft.icon')

# As of 0.22.0, served automatically from /assets/


###################################################################
if __name__ == '__main__':
  app.run_server(debug = True)
  # app.run_server()