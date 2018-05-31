###############################################
# Functions and variables for URL shortening
###############################################

chars = None
dna_to_code = dict()
code_to_dna = dict()

def __init_chars():
  global chars
  chars = [chr(s) for s in range(48, 48 + 10)] + [chr(s) for s in range(65, 65 + 26)] + [chr(s) for s in range(97, 97 + 26)]
  chars += ['-', '_', '~', '.']
  chars.remove('A')
  chars.remove('C')
  chars.remove('G')
  chars.remove('T')
  chars.remove('.')
  return

def __init_mappers():
  output = chars
  # All 3-mers of 61-length safe html character alphabet
  for idx in range(3-1):
    output = __append_alphabet(output, chars)
  triplets = output

  # All 8-mers DNA
  output = list('ACGT')
  for idx in range(8-1):
    output = __append_alphabet(output, list('ACGT'))
  eightmers = output

  global dna_to_code
  global code_to_dna
  for eightmer, triplet in zip(eightmers, triplets):
    dna_to_code[eightmer] = triplet
    code_to_dna[triplet] = eightmer
  return

def __append_alphabet(output, alphabet):
  new_output = []
  for o in output:
    for a in alphabet:
      new_output.append(o + a)
  return new_output

def parse_valid_url_path(url_path):
  ## Expected format:
  # [code][ACGT*].[indexnumber of cutsite]
  url_path = url_path.strip('/')
  if len(url_path) == 0:
    return False, None, None

  escape_chars = [s for s in list('ACGT.') if s in url_path]
  idx = min([url_path.index(s) for s in escape_chars])
  coded = url_path[:idx]
  tail = url_path[idx:]

  if len(coded) % 3 != 0:
    return False, None, None  

  seq = ''
  for jdx in range(0, len(coded), 3):
    w = coded[jdx : jdx + 3]
    seq += code_to_dna[w]


  if '.' not in tail:
    idx = len(tail)
  else:
    idx = tail.index('.')
  for ch in set(tail[:idx]):
    if ch not in list('ACGT'):
      return False, None, None
  seq += tail[:idx]

  if '.' not in tail:
    return True, seq, None
  try:
    cutsite_index = int(tail[idx+1:])
  except:
    return False, None, None
  return True, seq, cutsite_index

def encode_dna_to_url_path(seq, cutsite):
  code = ''
  for idx in range(0, len(seq), 8):
    chomp = seq[idx : idx + 8]
    if len(chomp) != 8:
      code += seq[idx:]
    else:
      code += dna_to_code[chomp]
  return '%s.%s' % (code, cutsite)

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
