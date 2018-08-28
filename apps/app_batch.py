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
        id = 'B_hidden-pred-df-stats',
        children = 'init'
      ),
      html.Div(
        id = 'B_hidden-cache-submit-button',
        children = '%s' % (time.time())
      ),
      html.Div(
        id = 'B_hidden-sort-module-interaction',
        children = '%s' % (time.time())
      ),
      html.Div(
        id = 'B_hidden-clickData',
        children = '%s init' % (time.time())
      ),
      html.Div(
        id = 'B_hidden-selected-id',
        children = ''
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

  ###################################################
  # Header
  ###################################################
  html.Div(
    [
      ###################################################
      # Upper header
      ###################################################
      header.get_navigation_header('batch'),

      ###################################################
      # Sequence box
      ###################################################
      html.Div([
        dcc.Textarea(
          id = 'B_textarea', 
          value = 'GCAATATCGTAGTCCGTCAAATTCAGCCCTGTTATCCCCGGCGTTATGTGTCAAATGGCGTAGAACTGGATTGACTGTTTGACGGTACCTGCTGATCGGTACGGTGACCGAGAATCTGTCGGGCTATGTCACTAATACTTT',
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
          html.Span(
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
                'Cutsite assumed 3nt upstream of PAM match. Supports IUPAC DNA encoding, ex: NNNRRT, NGG.',
                className = 'tooltiptext'
              ),
            ], 
            className = 'tooltip',
          ),
        ],
        style = dict(
          textAlign = 'center',
          marginBottom = 10,
        ),
      ),

      ###################################################
      # Advanced options: 
      ###################################################
      html.Div(
        [
          # header
          html.Div([
            html.Div(
              html.Strong(
                '▶ Advanced options',
                id = 'B_advanced_options_header_text',
              ),
              className = 'module_header_text'),
            ],
            id = 'B_advanced_options_header',
            style = dict(
              backgroundColor = 'rgba(0, 0, 0, 0.05)',
              height = 34,
              verticalAlign = 'middle',
            ),
          ),

          # Body
          html.Div(
            [
              html.Div(
                html.Strong('Evaluate gRNAs by...'),
                style = dict(
                  height = '50px',
                  lineHeight = '50px',
                  textAlign = 'left',
                  transform = 'translateX(35px)',
                ),
              ),

              # Row: Match sequence specification
              html.Div([
                html.Span(
                  'Repair frequency to a specific genotype',
                  style = dict(
                    textAlign = 'right',
                    lineHeight = 1.2,
                  ),
                  className = 'three columns',
                ),
                # Match sequence specification
                html.Div([
                  dcc.Textarea(
                    id = 'B_adv_matchseq', 
                    placeholder = 'Provide a DNA sequence',
                    style = dict(
                      fontFamily = 'monospace',
                      fontSize = 16,
                      resize = 'none',
                      height = '60px',
                      width = '95%',
                    ),
                  )],
                  style = dict(
                    verticalAlign = 'center',
                    whiteSpace = 'nowrap',
                  ),
                  className = 'nine columns'
                ),
                ##
                ], 
                style = dict(
                  marginBottom = 10,
                ),
                className = 'row',
              ),

              # Row: position of interest
              html.Div([
                html.Span(
                  'Distance to specific nucleotide',
                  style = dict(
                    textAlign = 'right',
                    lineHeight = 1.2,
                  ),
                  className = 'three columns',
                ),
                # Position of interest
                html.Div(
                  [
                    dcc.Input(
                      id = 'B_adv_position_of_interest',
                      type = 'text',
                      inputmode = 'numeric',
                      placeholder = '#',
                      min = 1,
                      step = 1,
                      style = dict(
                        width = 60,
                      ),
                    ),
                    html.Span(
                      id = 'B_adv_poi_selected_seq',
                    ),
                  ],
                  style = dict(
                  ),
                  className = 'nine columns',
                ),
                ##
                ],
                style = dict(
                  marginBottom = 15,
                ),
                className = 'row',
              ),

              # Row: deletion specification
              html.Div([
                html.Span(
                  'Frequency of deletions involving nucleotides',
                  style = dict(
                    textAlign = 'right',
                    lineHeight = 1.2,
                  ),
                  className = 'three columns',
                ),
                # Deletion specification
                html.Div(
                  [
                    dcc.Input(
                      id = 'B_adv_delstart',
                      type = 'text',
                      inputmode = 'numeric',
                      placeholder = '#',
                      min = 1,
                      step = 1,
                      style = dict(
                        width = 60,
                      ),
                    ),
                    html.Strong(
                      ' — '
                    ),
                    dcc.Input(
                      id = 'B_adv_delend',
                      type = 'text',
                      inputmode = 'numeric',
                      placeholder = '#',
                      min = 1,
                      step = 1,
                      style = dict(
                        width = 60,
                      ),
                    ),
                    html.Span(
                      id = 'B_adv_delseq',
                    ),
                  ],
                  style = dict(
                  ),
                  className = 'nine columns',
                ),
                ##
                ],
                style = dict(
                  marginBottom = 20,
                ),
                className = 'row',
              ),

              # Empty spacer
              html.Div(
                '',
                style = dict(
                  height = '10px',
                ),
              )

            ],
            id = 'B_advanced_options_body',
            style = dict(
              display = 'none',
            ),
            className = 'animate-top',
          ),
        ],
        id = 'B_advanced_options_module',
        style = dict(
          width = 750,
          margin = '0 auto',
          boxShadow = '1px 3px 6px 0 rgba(0, 0, 0, 0.2)',
          marginBottom = 30,
        )
      ),


      ###################################################
      # Click to run button + time estimate
      ###################################################
      # Time estimate
      html.P(
        id = 'B_estimated_runtime',
        children = 'Provide a sequence and PAM.',
        style = dict(
          textAlign = 'center',
        ),
      ),

      # Submit button
      html.Div([
        html.Button(
          'PREDICT REPAIR',
          id = 'B_submit_button',
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
              id = 'B_postcomp_module_header',
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
                  value = ['Cutsite', 'Precision', 'Frameshift (%)', 'Log phi', 'M.F. gt (%)'],
                  className = 'nine columns',
                ),
              ],
              style = dict(
                # width = '1050px',
                marginBottom = '5px',
                marginTop = '10px',
              ),
              className = 'row',
              id = 'B_row_dropdown-columns',
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
                  id = 'B_dropdown-sortcol',
                  options = [],
                  searchable = False,
                  clearable = False,
                  className = 'three columns',
                ),
                # Sort direction
                dcc.RadioItems(
                  id = 'B_sortdirection',
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
              id = 'B_row_dropdown-sortcol',
            ),

            # Links
            html.Div([
              html.Div(
                # Sharable link
                html.A(
                  '🔗 Shareable link to page before computation', 
                  id = 'B_page-link'
                )
              ),
              html.Div(
                # Download link: summary statistics
                html.A(
                  '📑 Download table of predictions',
                  id = 'B_download-link'
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
          id = 'B_hist-stats',
          config = dict(
            modeBarButtonsToRemove = modebarbuttons_2d,
            displaylogo = False,
            displayModeBar = False,
          ),
        ),
        id = 'B_hist-stats-div',
        style = dict(
          display = 'none',
          position = 'relative',
          zIndex = 1,
        )
      ),
    ],
    # body style
    id = 'B_postcomputation_settings',
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
          id = 'B_plot-stats',
          config = dict(
            modeBarButtonsToRemove = modebarbuttons_2d,
            displaylogo = False,
            displayModeBar = False,
          ),
        ),
        id = 'B_plot-stats-div',
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
  Output('B_hidden-cache-submit-button', 'children'),
  [Input('B_submit_button', 'n_clicks')])
def update_submit_button_time(n_clicks):
  return '%s' % (time.time())

@app.callback(
  Output('B_hidden-sort-module-interaction', 'children'),
  [Input('B_row_dropdown-columns', 'n_clicks'),
   Input('B_row_dropdown-sortcol', 'n_clicks')])
def update_sort_time(v1, v2):
  return '%s' % (time.time())

@app.callback(
  Output('B_hidden-clickData', 'children'),
  [Input('B_plot-stats', 'clickData')])
def update_hidden_clickdata(clickData):
  return '%s %s' % (time.time(), clickData['points'][0]['pointNumber'])

##
# URL and header callbacks
##
@app.callback(
  Output('B_textarea', 'value'),
  [Input('B_url', 'pathname')],
  [State('B_textarea', 'value')])
def update_textarea_from_url(url, default_value):
  valid_flag, dd = lib.parse_valid_url_path_batch(url)
  if valid_flag:
    return dd['seq']
  return default_value

@app.callback(
  Output('B_textbox_pam', 'value'),
  [Input('B_url', 'pathname')],
  [State('B_textbox_pam', 'value')])
def update_pam_from_url(url, default_value):
  valid_flag, dd = lib.parse_valid_url_path_batch(url)
  if valid_flag:
    return dd['pam']
  return default_value

@app.callback(
  Output('B_adv_matchseq', 'value'),
  [Input('B_url', 'pathname')],
  [State('B_adv_matchseq', 'value')])
def update_adv_matchseq_from_url(url, default_value):
  valid_flag, dd = lib.parse_valid_url_path_batch(url)
  if valid_flag:
    if dd['adv_seq_spec'] == '-':
      return default_value
    else:
      return dd['adv_seq_spec']
  return default_value

@app.callback(
  Output('B_adv_position_of_interest', 'value'),
  [Input('B_url', 'pathname')],
  [State('B_adv_position_of_interest', 'value')])
def update_adv_poi_from_url(url, default_value):
  valid_flag, dd = lib.parse_valid_url_path_batch(url)
  if valid_flag:
    if dd['adv_poi'] == '-':
      return default_value
    else:
      return dd['adv_poi']
  return default_value

@app.callback(
  Output('B_adv_delstart', 'value'),
  [Input('B_url', 'pathname')],
  [State('B_adv_delstart', 'value')])
def update_adv_delstart_from_url(url, default_value):
  valid_flag, dd = lib.parse_valid_url_path_batch(url)
  if valid_flag:
    if dd['adv_delstart'] == '-':
      return default_value
    else:
      return dd['adv_delstart']
  return default_value

@app.callback(
  Output('B_adv_delend', 'value'),
  [Input('B_url', 'pathname')],
  [State('B_adv_delend', 'value')])
def update_adv_delend_from_url(url, default_value):
  valid_flag, dd = lib.parse_valid_url_path_batch(url)
  if valid_flag:
    if dd['adv_delend'] == '-':
      return default_value
    else:
      return dd['adv_delend']
  return default_value

##
# Precomputation text / Advanced options callbacks
##
@app.callback(
  Output('B_estimated_runtime', 'children'),
  [Input('B_textarea', 'value'),
   Input('B_textbox_pam', 'value')])
def update_estimated_runtime(seq, pam):
  # Error catching
  if len(seq) < 70:
    return 'Error: Provide a sequence longer than 70 bp.'
  if len(seq) > 5000:
    return 'Error: Provide a sequence shorter than 5kb.'
  if len(pam) < 2 or len(pam) > 6:
    return 'Error: Provide a PAM between 2 and 6 bp long.'
  allowed_seq_chars = set(list('ACGTacgt'))
  for char in set(seq):
    if char not in allowed_seq_chars:
      return 'Error: Sanitize your sequence: %s disallowed' % (char)
  allowed_pam_chars = set(list('ACGTYRWSKMDVHBNacgtyrwskmdvhbn'))
  for char in set(pam):
    if char not in allowed_pam_chars:
      return 'Error: Sanitize your PAM: %s disallowed' % (char)
  if pam.count('N') == len(pam):
    return 'Error: PAM cannot only consist of N'


  pam_freq = lib.estimate_pam_freq(pam) * 2 # rc also
  num_est_pams = pam_freq * len(seq)
  est_time_per_pam = 0.2 # seconds
  est_runtime = est_time_per_pam * num_est_pams

  if est_runtime < 2:
    # ans = '1 second'
    crispr_words = ['CRISPR', 'gene editing', 'DNA', 'gene drive', 'non-homologous end-joining', 'homology-directed repair', 'microhomology-mediated end-joining', 'inDelphi', 'DNA microhomology', 'programmable nuclease', 'TAL effector nuclease', 'zinc-finger nuclease', 'genome editing', 'protospacer', 'protospacer-adjacent motif', 'prokaryotic antiviral defense mechanism', 'Streptococcus pyogenes', 'Cas9', 'single guide RNA', 'tracrRNA', 'crRNA', 'R loop', 'genetic engineering', 'gene knockout', 'computational biology', 'synthetic biology', 'disease correction', 'double-strand break']
    import random
    ans = 'faster than you can say "%s"' % (random.choice(crispr_words))
  elif est_runtime < 10:
    ans = '%s seconds' % (int(est_runtime))
  elif est_runtime < 60:
    ans = '%s0 seconds' % (int(round(est_runtime / 10)))
    if ans == '60 seconds':
      ans = '1 minute'
  elif est_runtime < 90:
    ans = '1 minute'
  elif est_runtime < 60*60:
    ans = '%s minutes' % (int(round(est_runtime / 60)))
  elif est_runtime < 1.5 * 60*60:
    ans = '1 hour'
  else:
    ans = '%s hours' % (int(round(est_runtime / (60*60))))
  if est_runtime > 25:
    # Address Heroku's 30-second timeout
    ans += '. Warning: Jobs over 30 seconds will time out.'
  return 'Estimated runtime: %s' % (ans)

@app.callback(
  Output('B_adv_poi_selected_seq', 'children'),
  [Input('B_adv_position_of_interest', 'value'),
   Input('B_textarea', 'value')])
def update_position_of_interest_selected_seq(poi, seq):
  # poi is 1-indexed
  poi_0idx = int(poi) - 1
  buff = 7
  if poi_0idx < buff or poi_0idx > len(seq) - buff:
    return ''
  selected_base = seq[poi_0idx]
  left = seq[poi_0idx - buff : poi_0idx]
  right = seq[poi_0idx + 1: poi_0idx + 1 + buff]
  def get_style_dict(color_char):
    return dict(
      fontFamily = 'monospace',
      fontSize = 14,
      color = '#%s' % (color_char * 6),
    )
  children = []
  gradient = list('EDCBA98')
  for nt, cc in zip(left, gradient):
    children.append(
      html.Span(nt, style = get_style_dict(cc)),
    )
  children.append(
    html.Strong(selected_base, style = get_style_dict('4')),
  )
  for nt, cc in zip(right, gradient[::-1]):
    children.append(
      html.Span(nt, style = get_style_dict(cc)),
    )
  return children

@app.callback(
  Output('B_adv_delseq', 'children'),
  [Input('B_adv_delstart', 'value'),
   Input('B_adv_delend', 'value'),
   Input('B_textarea', 'value')])
def update_selected_delseq(del_start, del_end, seq):
  # poi is 1-indexed, convert to 0-idx
  del_start = int(del_start) - 1
  del_end = int(del_end) - 1
  buff = 7
  if del_start >= del_end:
    return ''
  if del_start < buff or del_end > len(seq) - buff:
    return ''
  left = seq[del_start - buff : del_start]
  selected_bases = seq[del_start : del_end]
  right = seq[del_end : del_end + buff]
  def get_style_dict(color_char):
    return dict(
      fontFamily = 'monospace',
      fontSize = 14,
      color = '#%s' % (color_char * 6),
    )
  children = []
  gradient = list('EDCBA98')
  for nt, cc in zip(left, gradient):
    children.append(
      html.Span(nt, style = get_style_dict(cc)),
    )
  children.append(
    html.Strong(selected_bases, style = get_style_dict('4')),
  )
  for nt, cc in zip(right, gradient[::-1]):
    children.append(
      html.Span(nt, style = get_style_dict(cc)),
    )
  return children

##
# Submit button
##
@app.callback(
  Output('B_submit_button', 'children'),
  [Input('B_textarea', 'value'),
   Input('B_textbox_pam', 'value'),
   Input('B_estimated_runtime', 'children')])
def update_submit_button_text(seq, pam, est_runtime_text):
  if 'Error' in est_runtime_text:
    return 'PREDICT REPAIR'

  seq, pam = seq.upper(), pam.upper()
  num_grnas = 0
  seqs = [seq, lib.revcomp(seq)]
  cutsites = range(30, len(seq) - 30)
  for local_seq, grna_orient in zip(seqs, ['+', '-']):
    for local_cutsite in cutsites:
      cand_pam = local_seq[local_cutsite + 3 : local_cutsite + 3 + len(pam)]
      if lib.match(pam, cand_pam):
        num_grnas += 1
  return 'PREDICT REPAIR FOR %s gRNAs' % (num_grnas)

@app.callback(
  Output('B_submit_button', 'style'),
  [Input('B_estimated_runtime', 'children')],
  [State('B_submit_button', 'style')])
def update_submit_button_style(est_runtime_text, style):
  if 'Error' in est_runtime_text:
    style['backgroundColor'] = '#86898C'
    style['color'] = 'white'
  else:
    style['backgroundColor'] = '#00A0DC'
    style['color'] = 'white'
  return style


##
# Prediction callback
##
@app.callback(
  Output('B_hidden-pred-df-stats', 'children'),
  [Input('B_submit_button', 'n_clicks')],
  [State('B_textarea', 'value'),
   State('B_textbox_pam', 'value'),
   State('B_adv_matchseq', 'value'),
   State('B_adv_position_of_interest', 'value'),
   State('B_adv_delstart', 'value'),
   State('B_adv_delend', 'value'),
  ])
def update_pred_df_stats(nclicks, seq, pam, adv_matchseq, adv_poi, adv_delstart, adv_delend):
  # When submit button clicked, find all gRNAs matching PAM in sequence.
  # Advanced options:
  #   if matchseq is provided, include a column on
  #     sum frequencies of repair gts matching sequence
  #     e.g., pathogenic -> wildtype repair
  #   if deletion range is provided, include a column on
  #     sum frequencies of repair gts deleting specified positions.
  #   if position of interest is provided, include a column on
  #     cutsite distance to position of interest
  dd = defaultdict(list)
  all_stats = pd.DataFrame()

  assert pam.count('N') != len(pam)
  assert 2 <= len(pam) <= 6
  seq = seq.upper()
  pam = pam.upper()

  # Check and initialize advanced settings
  adv_matchseq_flag = False
  if adv_matchseq is not None and len(adv_matchseq) != 0:
    adv_matchseq = adv_matchseq.upper()
    adv_matchseq_flag = True
  adv_poi_flag = False
  if adv_poi is not None and len(adv_poi) > 0:
    # adv_poi is 1-indexed, switch to 0-index
    adv_poi = int(adv_poi) - 1
    adv_poi_flag = True
  adv_del_flag = False
  if adv_delstart is not None and adv_delend is not None:
    if len(adv_delstart) > 0 and len(adv_delend) > 0:
      adv_delstart, adv_delend = int(adv_delstart), int(adv_delend)
      if adv_delstart < adv_delend:
        adv_delstart -= 1
        adv_delend -= 1
        adv_del_flag = True


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
          cutsite_plus = local_cutsite
        else:
          cutsite_plus = len(seq) - local_cutsite
        dd['Cutsite'].append(cutsite_plus)

        # inDelphi predictions and standard statistics
        pred_df, stats = inDelphi.predict(local_seq, local_cutsite)
        all_stats = all_stats.append(stats, ignore_index = True)
        
        # Detailed link
        sm_link = lib.encode_dna_to_url_path_single(local_seq, local_cutsite)
        dd['URL'].append('https://dev.crisprindelphi.design%s' % (sm_link))

        # Handle advanced options
        if adv_matchseq_flag:
          inDelphi.add_genotype_column(pred_df, stats)
          crit = (pred_df['Genotype'] == adv_matchseq)
          matched_seq_freq = sum(pred_df[crit]['Predicted frequency'])
          dd['Repairs to spec.'].append(matched_seq_freq)

        if adv_poi_flag:
          if adv_poi > cutsite_plus:
            dist = abs(cutsite_plus - 1 - adv_poi)
          else:
            dist = abs(cutsite_plus - adv_poi)
          dd['Dist. to POI'].append(dist)

        if adv_del_flag:
          crit = (pred_df['Category'] == 'del') & (pred_df['Genotype position'] != 'e')
          delseq_freq = 0
          if grna_orient == '+':
            adv_delstart_local = adv_delstart
            adv_delend_local = adv_delend
          else:
            adv_delstart_local = len(seq) - adv_delend
            adv_delend_local = len(seq) - adv_delstart
          for jdx, row in pred_df[crit].iterrows():
            del_start = local_cutsite - row['Length'] + row['Genotype position']
            del_end = del_start + row['Length']
            if del_start <= adv_delstart_local < adv_delend_local <= del_end:
              delseq_freq += row['Predicted frequency']
          dd['Deletes spec.'].append(delseq_freq)

  # Add metadata columns and advanced settings
  for col in dd:
    all_stats[col] = dd[col]

  # Switch phi to log phi
  all_stats['Log phi'] = np.log(all_stats['Phi'])  
  all_stats = all_stats.drop(['Phi'], axis = 1)

  # Sort by cutsite and relabel indices
  all_stats = all_stats.sort_values(by = 'Cutsite')
  all_stats = all_stats.reset_index(drop = True)

  all_stats['ID'] = all_stats.index + 1

  return all_stats.to_csv()
  # return (pred_df.to_csv(), pd.DataFrame(stats, index = [0]).to_csv())

##
# Module header callbacks, Advanced options hiding/showing
##
@app.callback(
  Output('B_postcomp_module_header', 'children'),
  [Input('B_hidden-pred-df-stats', 'children')],
  [State('B_textarea', 'value'),
   State('B_textbox_pam', 'value')])
def update_postcomp_module_header(all_stats_string, seq, pam):
  if all_stats_string == 'init':
    assert False, 'init'
  stats = pd.read_csv(StringIO(all_stats_string), index_col = 0)
  return 'Results of %s gRNAs with %s PAM found in %s-bp query' % (len(stats), pam, len(seq))

@app.callback(
  Output('B_advanced_options_body', 'style'),
  [Input('B_advanced_options_header', 'n_clicks'),
   Input('B_url', 'pathname')],
  [State('B_advanced_options_body', 'style')])
def update_adv_options_body_style(n_clicks, url, prev_style):
  new_style = prev_style
  if n_clicks is None:
    valid_flag, dd = lib.parse_valid_url_path_batch(url)
    if valid_flag and dd['adv_flag'] == True:
      del new_style['display']

  elif n_clicks > 0:  # ignore first automatic click triggered by page load
    if 'display' in prev_style:
      del new_style['display']
    else:
      new_style['display'] = 'none'
  return new_style

@app.callback(
  Output('B_advanced_options_header_text', 'children'),
  [Input('B_advanced_options_header', 'n_clicks')],
  [State('B_advanced_options_header_text', 'children')])
def update_adv_options_header_text(n_clicks, prev_text):
  if n_clicks > 0:
    if '▶' in prev_text:
      new_arrow = '▼'
    else:
      new_arrow = '▶'
  return '%s Advanced options' % (new_arrow)

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

@app.callback(
  Output('B_dropdown-sortcol', 'value'),
  [Input('B_dropdown-sortcol', 'options')],
  [State('B_url', 'pathname'),
   State('B_dropdown-sortcol', 'value'),
   State('B_advanced_options_module', 'n_clicks'),
   State('B_row_dropdown-columns', 'n_clicks'),
   State('B_row_dropdown-sortcol', 'n_clicks'),
   ])
def update_sortcol_value_from_url(options, url, prev_value, nc1, nc2, nc3):
  if nc1 or nc2 or nc3:
    # If clicked on any module that might change the sortcol
    return prev_value
  valid_flag, dd = lib.parse_valid_url_path_batch(url)
  if not valid_flag or dd['sort_by'] == '-':
    return prev_value
  else:
    all_options = [s['value'] for s in options]
    idx = int(dd['sort_by'])
    return sorted(all_options)[idx]

@app.callback(
  Output('B_dropdown-columns', 'options'),
  [Input('B_hidden-pred-df-stats', 'children')],
  [State('B_dropdown-columns', 'options')]
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
  Output('B_dropdown-columns', 'value'),
  [Input('B_dropdown-columns', 'options')],
  [State('B_dropdown-columns', 'value'),
   State('B_url', 'pathname'),
   State('B_row_dropdown-columns', 'n_clicks')]
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

@app.callback(
  Output('B_sortdirection', 'value'),
  [Input('B_dropdown-sortcol', 'options')],
  [State('B_url', 'pathname'),
   State('B_sortdirection', 'value')])
def update_sortdir_from_url(sort_options, url, prev_value):
  valid_flag, dd = lib.parse_valid_url_path_batch(url)
  if valid_flag:
    return dd['sort_dir']
  else:
    return prev_value

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
  if all_stats_string == 'init':
    assert False, 'init'
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
    if sort_direction == 'Ascending':
      ascending_flag = True
    else:
      ascending_flag = False
    stats = stats.sort_values(by = sort_col, ascending = ascending_flag)

  # Reformat floats
  stats_cols = list(stats.columns)
  nonstat_cols = ['ID', 'gRNA', 'gRNA orientation', 'PAM', 'URL']
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
  [Input('B_hidden-clickData', 'children'),
   Input('B_hidden-cache-submit-button', 'children'),
   Input('B_dropdown-columns', 'value'),
   Input('B_dropdown-sortcol', 'value'),
   Input('B_table-stats', 'rows')],
  [State('B_table-stats', 'selected_row_indices'),
   State('B_hidden-sort-module-interaction', 'children'),
   State('B_hidden-selected-id', 'children'),
   State('B_url', 'pathname'),
   State('B_postcomputation_settings', 'n_clicks'),
   State('B_plot-stats-div', 'n_clicks'),
   State('B_submit_button', 'n_clicks'),
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
    new_idx = int(df[df['ID'] == int(prev_id)].index[0])
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
  Output('B_hidden-selected-id', 'children'),
  [Input('B_table-stats', 'selected_row_indices')],
  [State('B_table-stats', 'rows')])
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
  Output('B_plot-stats-div', 'style'),
  [Input('B_plot-stats', 'figure')])
def update_stats_plot_style(fig):
  if fig is None:
    return {'display': 'none'}
  else:
    return {}

@app.callback(
  Output('B_hist-stats-div', 'style'),
  [Input('B_hist-stats', 'figure')])
def update_hist_plot_style(fig):
  if fig is None:
    return {'display': 'none'}
  else:
    return {}

@app.callback(
  Output('B_postcomputation_settings', 'style'),
  [Input('B_plot-stats', 'figure')])
def update_postcomputation_settings_style(fig):
  if fig is None:
    return {'display': 'none'}
  else:
    return {}

########################################################
# Plot stats callback
########################################################
@app.callback(
    Output('B_plot-stats', 'figure'),
    [Input('B_table-stats', 'rows'),
     Input('B_table-stats', 'selected_row_indices')])
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

  # Format y tick texts: ID, gRNA, PAM, orientation, URL.
  yticktexts = []
  fixedwidth_ids = lib.get_fixedwidth_ID(df['ID'])
  for idx, row in df.iterrows():
    row_text = '%s %s %s <a href="%s">details</a> %s' % (row['gRNA'], row['PAM'], row['gRNA orientation'], row['URL'], fixedwidth_ids[idx])
    yticktexts.append(row_text)


  # Subplot formatting
  fig['layout']['yaxis1'].update(
    fixedrange = True,
    tickvals = np.arange(1, len(df.index) + 1),
    ticktext = yticktexts,
    tickfont = dict(
      size = 12,
      family = 'monospace',
    ),
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
    Output('B_hist-stats', 'figure'),
    [Input('B_table-stats', 'rows'),
     Input('B_table-stats', 'selected_row_indices')])
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


##
# Download callbacks
##
@app.callback(
  Output('B_download-link', 'href'), 
  [Input('B_table-stats', 'rows')])
def update_link(rows):
  df = pd.DataFrame(rows)
  stats_cols = list(df.columns)
  nonstat_cols = ['gRNA', 'gRNA orientation', 'PAM', 'URL', 'ID']
  for nonstat_col in nonstat_cols:
    stats_cols.remove(nonstat_col)
  df = df[nonstat_cols + lib.order_chosen_columns(stats_cols)]

  time = str(datetime.datetime.now()).replace(' ', '_').replace(':', '-')
  link_fn = '/dash/urlToDownloadBatch?value={}'.format(time)
  df.to_csv('user-csvs/%s.csv' % (time), index = False)
  return link_fn

##
# Flask serving
##
@app.server.route('/dash/urlToDownloadBatch') 
def download_csv_batch():
  value = flask.request.args.get('value')
  # create a dynamic csv or file here using `StringIO` 
  # (instead of writing to the file system)
  local_csv_fn = value.split('/')[-1]
  return flask.send_file(
    open('user-csvs/%s.csv' % (local_csv_fn), 'rb'),
    mimetype = 'text/csv',
    attachment_filename = 'inDelphiBatch_output.csv',
    as_attachment = True,
  )


##
# Page link callback
##
@app.callback(
  Output('B_page-link', 'href'),
  [Input('B_textarea', 'value'),
   Input('B_textbox_pam', 'value'),
   Input('B_advanced_options_body', 'style'),
   Input('B_adv_matchseq', 'value'),
   Input('B_adv_position_of_interest', 'value'),
   Input('B_adv_delstart', 'value'),
   Input('B_adv_delend', 'value'),
   Input('B_dropdown-columns', 'value'),
   Input('B_dropdown-columns', 'options'),
   Input('B_dropdown-sortcol', 'value'),
   Input('B_sortdirection', 'value'),
   Input('B_table-stats', 'selected_row_indices'),
  ])
def update_pagelink(textarea, pam, adv_style, adv_seq_spec, adv_poi, adv_delstart, adv_delend, chosen_columns, column_options, sort_by, sort_dir, selected_row):
  adv_flag = bool('display' not in adv_style)
  url = 'https://dev.crisprindelphi.design%s' % (lib.encode_dna_to_url_path_batch(textarea, pam, adv_flag, adv_seq_spec, adv_poi, adv_delstart, adv_delend, chosen_columns, column_options, sort_by, sort_dir, selected_row))
  return url