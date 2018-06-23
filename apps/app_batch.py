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

###################################################################
###################################################################
##
# App layout
##
layout = html.Div([
  html.Div([

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
                  'Cutsite assumed 3nt upstream of PAM match. Supports IUPAC DNA encoding, ex: NNNRRT, NGG.',
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
        # Advanced options: Match sequence
        ###################################################
        html.Div([
          dcc.Textarea(
            id = 'B_adv_matchseq', 
            placeholder = 'Provide a DNA sequence to be matched against all repair genotypes at all gRNAs.',
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

        # Position of interest
        html.Div(
          [
            dcc.Input(
              id = 'B_adv_position_of_interest',
              type = 'number',
              inputmode = 'numeric',
              min = 1,
              step = 1,
            ),
            html.Span(
              id = 'B_adv_poi_selected_seq',
            ),
          ],
        ),

        # positions to be deleted

        ###################################################
        # Click to run button + time estimate
        ###################################################
        html.P(
          id = 'B_estimated_runtime',
          children = 'Provide a sequence and PAM.',
          style = dict(
            textAlign = 'center',
          ),
        ),

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
        backgroundColor = 'white',
        width = '1010px',
      ),
    ),

    ###################################################
    # Post-computation settings
    ###################################################
    html.Div(
      [
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

        # Sharable link, move this later
        html.Div(
          html.A(
            '🔗 Shareable link to page before computation', 
            id = 'B_page-link'
          ),
        ),

        # Download link: summary statistics
        html.Div(
          html.A(
            '📑 Download table of predictions',
            id = 'B_download-link'
          ),
        ),

        # Hists
        html.Div(
          dcc.Graph(
            id = 'B_hist-stats',
            config = dict(
              modeBarButtonsToRemove = modebarbuttons_2d,
              displaylogo = False,
            ),
          ),
          id = 'B_hist-stats-div',
          style = dict(
            display = 'none',
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
          )
        ),

      ],
      # body style
      style = dict(
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
# Hidden button callbacks
##
@app.callback(
  Output('B_hidden-cache-submit-button', 'children'),
  [Input('B_submit_button', 'n_clicks')])
def update_submit_button_time(n_clicks):
  return '%s' % (time.time())

@app.callback(
  Output('B_hidden-sort-module-interaction', 'children'),
  [Input('B_sortdirection', 'value'),
   Input('B_dropdown-sortcol', 'value')])
def update_submit_button_time(v1, v2):
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
# Precomputation text / Advanced options callbacks
##
@app.callback(
  Output('B_estimated_runtime', 'children'),
  [Input('B_textarea', 'value'),
   Input('B_textbox_pam', 'value')])
def update_estimated_runtime(seq, pam):
  # Error catching
  if len(seq) < 70:
    return 'Provide a sequence longer than 70 bp.'
  if len(seq) > 5000:
    return 'Provide a sequence shorter than 5kb.'
  if len(pam) < 2 or len(pam) > 6:
    return 'Provide a PAM between 2 and 6 bp long.'
  allowed_seq_chars = set(list('ACGTacgt'))
  for char in set(seq):
    if char not in allowed_seq_chars:
      return 'Sanitize your sequence: %s disallowed' % (char)
  allowed_pam_chars = set(list('ACGTYRWSKMDVHBNacgtyrwskmdvhbn'))
  for char in set(pam):
    if char not in allowed_pam_chars:
      return 'Sanitize your PAM: %s disallowed' % (char)
  if pam.count('N') == len(pam):
    return 'PAM cannot only consist of N'


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
  return 'Estimated runtime: %s' % (ans)

@app.callback(
  Output('B_adv_poi_selected_seq', 'children'),
  [Input('B_adv_position_of_interest', 'value'),
   Input('B_textarea', 'value')])
def update_position_of_interest_selected_seq(poi, seq):
  # poi is 1-indexed
  poi_0idx = poi - 1
  buff = 8
  if poi_0idx < buff or poi_0idx > len(seq) - buff:
    return ''
  selected_base = seq[poi_0idx]
  left = seq[poi_0idx - buff : poi_0idx - 1]
  right = seq[poi_0idx + 1: poi_0idx + buff]
  def get_style_dict(color_char):
    return dict(
      fontFamily = 'monospace',
      fontSize = 14,
      color = '#%s' % (color_char * 6),
    )
  children = []

  gradient = list('EDCBA987')
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
  ])
def update_pred_df_stats(nclicks, seq, pam, adv_matchseq, adv_poi):
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
  if adv_poi is not None:
    adv_poi = int(adv_poi)
    adv_poi_flag = True

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

        pred_df, stats = inDelphi.predict(local_seq, local_cutsite)
        all_stats = all_stats.append(stats, ignore_index = True)
        
        sm_link = lib.encode_dna_to_url_path_single(local_seq, local_cutsite)
        dd['URL'].append('https://dev.crisprindelphi.design%s' % (sm_link))

        if adv_matchseq_flag:
          inDelphi.add_genotype_column(pred_df, stats)
          crit = (pred_df['Genotype'] == adv_matchseq)
          matched_seq_freq = sum(pred_df[crit]['Predicted frequency'])
          dd['Matched repair'].append(matched_seq_freq)

        if adv_poi_flag:
          # adv_poi is 1-indexed, switch to 0-index
          adv_poi_0idx = adv_poi - 1
          if adv_poi < cutsite_plus:
            dist = abs(cutsite_plus - 1 - adv_poi_0idx)
          else:
            dist = abs(cutsite_plus - adv_poi_0idx)
          dd['Dist. to POI'].append(dist)

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
  Output('B_dropdown-columns', 'options'),
  [Input('B_hidden-pred-df-stats', 'children')],
  [State('B_dropdown-columns', 'options')]
  )
def update_columns_options(all_stats_string, prev_options):
  stats = pd.read_csv(StringIO(all_stats_string), index_col = 0)
  options = prev_options

  for d in ['Matched repair', 'Dist. to POI']:
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
  [Input('B_hidden-pred-df-stats', 'children')],
  [State('B_dropdown-columns', 'value')]
  )
def update_columns_value(all_stats_string, prev_value):
  stats = pd.read_csv(StringIO(all_stats_string), index_col = 0)
  value = prev_value

  for td in ['Matched repair', 'Dist. to POI']:
    if td in stats.columns:
      if td not in value:
        value.append(td)
    else:
      if td in value:
        value.remove(td)
  return value

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
    if sort_direction == 'Up':
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
   Input('B_hidden-sort-module-interaction', 'children'),
   Input('B_table-stats', 'rows')],
  [State('B_table-stats', 'selected_row_indices'),
   State('B_hidden-selected-id', 'children')])
def update_statstable_selected(clickData, submit_time, sort_time, rows, selected_row_indices, prev_id):
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

  if sort_intxn:
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
  fig['layout']['showlegend'] = False
  fig['layout']['width'] = 275 + len(stats_cols) * 150
  fig['layout']['height'] = 140
  fig['layout']['margin'] = {
    'l': 250,
    'r': 25,
    't': 60,
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
  ])
def update_pagelink(textarea, pam):
  return 'https://dev.crisprindelphi.design%s' % (lib.encode_dna_to_url_path_batch(textarea, pam))