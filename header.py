import dash
import dash_core_components as dcc
import dash_html_components as html

divider_text = '   •   '

def get_navigation_header(page_nm):
  font_size_param = 16

  dot_style = dict(
    color = 'gray',
    fontSize = '%spx' % (font_size_param),
  )

  default_style = dict(
    position = 'relative',
    textDecoration = 'none',
    textTransform = 'uppercase',
    fontFamily = 'sans-serif',
    color = 'black',
    marginBottom = '2px',
    fontSize = '%spx' % (font_size_param),
  )

  import copy
  selected_style = copy.copy(default_style)
  selected_style['borderBottom'] = '1px solid'

  styles = dict()
  for nm in ['single', 'batch', 'gene', 'guide', 'about']:
    if page_nm == nm:
      styles[nm] = selected_style
    else:
      styles[nm] = default_style

  return html.Div(
    [
      html.H4(
        'inDelphi',
        style = dict(
          textAlign = 'center',
        ),
      ),

      html.Div(
        [
          html.A(
            'Single mode',
            href = 'https://www.crisprindelphi.design/single',
            style = styles['single'],
            className = 'dynamicunderline',
          ),
          html.Span(
            divider_text,
          ),
          html.A(
            'Batch mode',
            href = 'https://www.crisprindelphi.design/batch',
            style = styles['batch'],
            className = 'dynamicunderline',
          ),
          html.Span(
            divider_text,
          ),
          html.A(
            'Gene mode',
            href = 'https://www.crisprindelphi.design/gene',
            style = styles['gene'],
            className = 'dynamicunderline',
          ),
          html.Span(
            divider_text,
          ),
          html.A(
            'User guide',
            href = 'https://www.crisprindelphi.design/guide',
            style = styles['guide'],
            className = 'dynamicunderline',
          ),
          html.Span(
            divider_text,
          ),
          html.A(
            'About',
            href = 'https://www.crisprindelphi.design/about',
            style = styles['about'],
            className = 'dynamicunderline',
          ),
        ],
        style = dict(
          marginBottom = 20,
          textAlign = 'center',
        ),
        className = 'row',
      ),

    ],
  )