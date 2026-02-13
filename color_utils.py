import colorsys

# =====================
# カラー設定とグラデーション用関数
# =====================

def lighten_color(hex_color, factor):
    hex_color = hex_color.lstrip('#')
    r, g, b = tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))
    h, l, s = colorsys.rgb_to_hls(r / 255, g / 255, b / 255)
    l = min(1.0, l + factor * (1.0 - l))
    r_new, g_new, b_new = colorsys.hls_to_rgb(h, l, s)
    return '#{:02x}{:02x}{:02x}'.format(int(r_new * 255), int(g_new * 255), int(b_new * 255))

def get_or_create_gradient(dwg, base_color, grad_dict):
    if base_color in grad_dict:
        return grad_dict[base_color]

    grad_id = f'grad_{len(grad_dict)}'
    light_color = lighten_color(base_color, 0.4)
    lighter_color = lighten_color(base_color, 0.7)

    grad = dwg.linearGradient(start=('0%', '100%'), end=('0%', '0%'), id=grad_id)
    grad.add_stop_color(offset='0.0', color=base_color)
    grad.add_stop_color(offset='0.5', color=light_color)
    grad.add_stop_color(offset='1.0', color=lighter_color)
    dwg.defs.add(grad)

    grad_dict[base_color] = grad_id
    return grad_id

def get_domain_color(domain_name, color_map, palette):
    """
    domain名ごとに一貫した色を返す
    """
    if domain_name not in color_map:
        color_map[domain_name] = palette[len(color_map) % len(palette)]
    return color_map[domain_name]
