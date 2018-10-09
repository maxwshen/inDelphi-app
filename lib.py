import numpy as np

###############################################
# Functions and variables for URL shortening
###############################################

chars = None
dna_to_code = dict()
code_to_dna = dict()

KMER_LEN = 9

def __init_chars():
  global chars
  chars = [chr(s) for s in range(48, 48 + 10)] + [chr(s) for s in range(65, 65 + 26)] + [chr(s) for s in range(97, 97 + 26)]
  chars += ['-', '_', '~', '.']
  chars.remove('_')
  return

def __init_mappers():
  output = chars
  # All 3-mers of 65-length safe html character alphabet
  for idx in range(3-1):
    output = __append_alphabet(output, chars)
  triplets = output

  # All 9-mers DNA
  output = list('ACGT')
  for idx in range(KMER_LEN-1):
    output = __append_alphabet(output, list('ACGT'))
  kmers = output

  global dna_to_code
  global code_to_dna
  for kmer, triplet in zip(kmers, triplets):
    dna_to_code[kmer] = triplet
    code_to_dna[triplet] = kmer
  return

def __append_alphabet(output, alphabet):
  new_output = []
  for o in output:
    for a in alphabet:
      new_output.append(o + a)
  return new_output

def parse_coded_seq_leftover(dd, coded_nm, leftover_nm):
  # Process encoded DNA
  if len(dd[coded_nm]) != 1 and len(dd[coded_nm]) % 3 != 0:
    return '-'
  if dd[coded_nm] == '-':
    return dd[leftover_nm]

  seq = ''
  for jdx in range(0, len(dd[coded_nm]), 3):
    w = dd[coded_nm][jdx : jdx + 3]
    seq += code_to_dna[w]
  if dd[leftover_nm] != '-':
    seq += dd[leftover_nm]
  return seq

def encode_dna(seq):
  if seq is None or len(seq) == 0:
    return '-', '-'

  if len(seq) < KMER_LEN:
    return '-', seq

  encodeddna = ''
  for idx in range(0, len(seq), KMER_LEN):
    chomp = seq[idx : idx + KMER_LEN]
    if len(chomp) == KMER_LEN:
      encodeddna += dna_to_code[chomp]
    else:
      break
  if len(seq[idx:]) != KMER_LEN:
    leftoverdna = seq[idx:]
  else:
    leftoverdna = '-'
  return encodeddna, leftoverdna

###############################################
# Single
###############################################

def parse_valid_url_path_single(url_path):
  ## Expected format:
  # [celltype]_[encodedDNA]_[leftoverDNA]_[cutsite]
  if url_path[:len('/single_')] != '/single_':
    return False, None, None, None

  url_path = url_path.replace('/single_', '')
  if len(url_path) == 0 or '_' not in url_path:
    return False, None, None, None

  parts = url_path.split('_')
  cats = ['celltype', 'coded', 'leftover', 'cutsite']
  if len(parts) != len(cats):
    return False, None, None, None
  dd = dict()
  for idx, cat in enumerate(cats):
    dd[cat] = parts[idx]

  seq = parse_coded_seq_leftover(dd, 'coded', 'leftover')
  return True, celltype, seq, int(dd['cutsite'])

def encode_dna_to_url_path_single(seq, cutsite, celltype):
  seq = seq.upper()
  encodeddna, leftoverdna = encode_dna(seq)
  return '/single_%s_%s_%s_%s' % (celltype, encodeddna, leftoverdna, cutsite)


###############################################
# Batch
###############################################

def parse_valid_url_path_batch(url_path):
  ## Expected format:
  # [encodedDNA]_[leftoverDNA]_[pam in plaintext] + more
  dd = dict()
  if url_path[:len('/batch_')] != '/batch_':
    return False, dd

  url_path = url_path.replace('/batch_', '')
  if len(url_path) == 0 or '_' not in url_path:
    return False, dd

  parts = url_path.split('_')
  cats = ['coded', 'leftover', 'pam', 'adv_flag', 'coded_spec', 'leftover_spec', 'adv_poi', 'adv_delstart', 'adv_delend', 'chosen_columns', 'sort_by', 'sort_dir', 'row_select']
  if len(parts) != len(cats):
    return False, dd
  for idx, cat in enumerate(cats):
    dd[cat] = parts[idx]

  dd['seq'] = parse_coded_seq_leftover(dd, 'coded', 'leftover')
  dd['adv_seq_spec'] = parse_coded_seq_leftover(dd, 'coded_spec', 'leftover_spec')

  # Reword some values
  if dd['adv_flag'] == '1':
    dd['adv_flag'] = True
  elif dd['adv_flag'] == '0':
    dd['adv_flag'] = False

  if dd['sort_dir'] == '1':
    dd['sort_dir'] = 'Ascending'
  else:
    dd['sort_dir'] = 'Descending'

  return True, dd

def encode_dna_to_url_path_batch(seq, pam, adv_flag, adv_seq_spec, adv_poi, adv_delstart, adv_delend, chosen_columns, column_options, sort_by, sort_dir, selected_row):
  seq, pam = seq.upper(), pam.upper()
  edna, ldna = encode_dna(seq)
  edna2, ldna2 = encode_dna(adv_seq_spec)

  if adv_flag == True:
    adv_flag_val = '1'
  else:
    adv_flag_val = '0'

  adv_poi = transform_empty_value_to_dash(adv_poi)
  adv_delstart = transform_empty_value_to_dash(adv_delstart)
  adv_delend = transform_empty_value_to_dash(adv_delend)
  sort_by = transform_empty_value_to_dash(sort_by)

  binary_flags_chosen_cols = ''
  for co in sorted([s['value'] for s in column_options]):
    if co in chosen_columns:
      binary_flags_chosen_cols += '1'
    else:
      binary_flags_chosen_cols += '0'

  if sort_by != '-':
    sort_by = sorted(chosen_columns).index(sort_by)

  if sort_dir == 'Ascending':
    sort_dir_val = '1'
  else:
    sort_dir_val = '0'

  if selected_row == []:
    selected_row_val = '-'
  else:
    selected_row_val = selected_row[0]

  items = [
    edna, 
    ldna, 
    pam, 
    adv_flag_val, 
    edna2, 
    ldna2, 
    adv_poi, 
    adv_delstart, 
    adv_delend, 
    binary_flags_chosen_cols, 
    sort_by, 
    sort_dir_val, 
    selected_row_val
  ]
  return '/batch_%s' % ('_'.join([str(s) for s in items]))

def transform_empty_value_to_dash(val):
  if val is None or len(val) == 0 or val == 'None':
    return '-'
  else:
    return val

__init_chars()
__init_mappers()

###############################################
# Compbio operations
###############################################

def revcomp(seq):
  rc_mapper = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'}
  rc_seq = []
  for c in seq:
    if c in rc_mapper:
      rc_seq.append(rc_mapper[c])
    else:
      rc_seq.append(c)
  return ''.join(rc_seq[::-1])

def pam_shift(text1, text2, text_pam, direction):
  seq = text1 + text2
  cutsite = len(text1)

  if direction == 'right':
    cutsites = range(cutsite + 1, len(seq))
  elif direction == 'left':
    cutsites = range(cutsite - 1, 0, -1)

  for ct in cutsites:
    candidate_pam = seq[ct + 3 : ct + 6]
    if match(text_pam, candidate_pam):
      return seq[:ct], seq[ct:]
  return None

mapper = {
  'A': list('A'),
  'C': list('C'),
  'G': list('G'),
  'T': list('T'),
  'Y': list('CT'),
  'R': list('AG'),
  'W': list('AT'),
  'S': list('GC'),
  'K': list('TG'),
  'M': list('AC'),
  'D': list('AGT'),
  'V': list('ACG'),
  'H': list('ACT'),
  'B': list('CGT'),
  'N': list('ACGT'),
}
def match(template, dna):
  if len(dna) != len(template):
    return False
  for char, t in zip(dna, template):
    if char not in mapper[t]:
      return False
  return True

def estimate_pam_freq(pam):
  factor = 1
  for char in pam:
    factor *= ( len(mapper[char]) / 4)
  return factor

###############################################
# Colors
###############################################

def get_color(stats_col):
  if stats_col in ['Cutsite', 'Dist. to POI']:
    return '#86898C'
  if stats_col == 'Exp. indel len':
    return '#86898C'
  if stats_col == 'Frame +0 (%)':
    return '#68C7EC'
  if stats_col == 'Frame +1 (%)':
    return '#68C7EC'
  if stats_col == 'Frame +2 (%)':
    return '#68C7EC'
  if stats_col == 'Frameshift (%)':
    return '#00A0DC'
  if stats_col == 'M.F. del (%)':
    return '#ED4795'
  if stats_col == 'M.F. ins (%)':
    return '#F47B16'
  if stats_col == 'M.F. gt (%)':
    return '#7CB82F'
  if stats_col == 'Log phi':
    return '#EC4339'
  if stats_col == 'Precision':
    return '#00AEB3'
  if stats_col in ['Repairs to spec.', 'Deletes spec.']:
    return '#C11F1D'
  return '#333333' # default

###############################################
# Batch mode: xaxis ticks 
###############################################
def get_batch_statcol_xrange(stats, stat_nm):
  if '(%)' in stat_nm:
    buff = 3
  elif stat_nm == 'Exp. indel len':
    buff = 1
  elif stat_nm == 'Log phi':
    buff = 0.1
  elif stat_nm == 'Precision':
    buff = 0.05 
  elif stat_nm == 'Cutsite':
    buff = 10
  elif stat_nm in ['Repairs to spec.', 'Deletes spec.']:
    buff = 5
  elif stat_nm == 'Dist. to POI':
    buff = 5
  else: # default
    buff = 0
  return [min(stats) - buff, max(stats) + buff]

# def get_batch_statcol_xticks(stats):
  # pass
  # return

def get_batch_select_line(x0 = 0, x1 = 0, y0 = 0, y1 = 0, xref = '', yref = ''):
  return dict(
    type = 'line',
    xref = xref,
    yref = yref,
    x0 = x0,
    x1 = x1,
    y0 = y0,
    y1 = y1,
    opacity = 0.8,
    line = dict(
      color = 'rgb(33, 33, 33)',
      width = 1,
      dash = 'dot',
    )
  )

def rename_batch_columns(stats):
  name_changes = {
    'Frameshift frequency': 'Frameshift (%)',
    'Frame +0 frequency': 'Frame +0 (%)',
    'Frame +1 frequency': 'Frame +1 (%)',
    'Frame +2 frequency': 'Frame +2 (%)',
    'Highest outcome frequency': 'M.F. gt (%)',
    'Highest del frequency': 'M.F. del (%)',
    'Highest ins frequency': 'M.F. ins (%)',
    'Expected indel length': 'Exp. indel len',
  }
  for col in stats:
    if col in name_changes:
      stats[name_changes[col]] = stats[col]
      stats.drop([col], axis = 1, inplace = True)
  return stats

def order_chosen_columns(cols):
  preferred_order = [
    'Cutsite',
    'Dist. to POI',
    'Repairs to spec.',
    'Deletes spec.',
    'Precision',
    'Frameshift (%)',
    'Frame +0 (%)',
    'Frame +1 (%)',
    'Frame +2 (%)',
    'Log phi',
    'M.F. gt (%)',
    'M.F. del (%)',
    'M.F. ins (%)',
    'Exp. indel len',
  ]
  reordered = []
  for pref in preferred_order:
    if pref in cols:
      reordered.append(pref)
  return reordered

def get_x_domains(num_cols):
  # Ensure uniform and consistent horizontal spacing with variable number of columns
  margin_pct = 0.12

  domains = []
  for leftside in np.arange(0, 1, 1/num_cols):
    size = 1 / num_cols
    margin_size = size * margin_pct 
    rightside = leftside + size
    domains.append([leftside + margin_size, rightside - margin_size])
  return domains

def get_fixedwidth_ID(ids):
  largest_len = len(str(max(ids)))
  fw_ids = []
  for item in ids:
    num_spaces = largest_len - len(str(item))
    fw_id = '%s#%s' % (' ' * num_spaces, item)
    fw_ids.append(fw_id)
  return fw_ids