#!/usr/bin/env python3
"""
Statistical significance utilities for plotting
"""

def format_significance_label(p_value):
    """Convert a p value bucket into the familiar star annotation."""
    if p_value is None:
        return ''
    if p_value < 0.001:
        return '***'
    if p_value < 0.01:
        return '**'
    if p_value < 0.05:
        return '*'
    return 'n.s.'


def draw_significance_annotations(ax, names, positions, means, sems, annotations,
                                   min_p_value=0.05):
    """Draw horizontal significance bars with stars above the requested bars."""
    if not annotations or not names:
        return

    name_to_pos = dict(zip(names, positions))
    name_to_top = {
        name: mean + sem for name, mean, sem in zip(names, means, sems)
    }
    lower, upper = ax.get_ylim()
    y_range = max(upper - lower, 1e-3)
    spacing = max(y_range * 0.08, 0.05)
    base_value = max(name_to_top.values(), default=upper)
    next_y = max(base_value + spacing, upper)
    max_y = next_y

    for idx, annotation in enumerate(annotations):
        group_names = annotation.get('groups') or annotation.get('pair') or annotation.get('pairs')
        if not group_names or len(group_names) != 2:
            continue
        name_a, name_b = group_names
        pos_a = name_to_pos.get(name_a)
        pos_b = name_to_pos.get(name_b)
        if pos_a is None or pos_b is None or pos_a == pos_b:
            continue

        x1, x2 = sorted((pos_a, pos_b))
        y = next_y
        p_value = annotation.get('p_value')
        if p_value is not None and p_value >= min_p_value and not annotation.get('force_draw'):
            continue

        text_label = annotation.get('text') or annotation.get('label')
        if not text_label:
            text_label = format_significance_label(p_value)
        if not text_label:
            continue

        tick_height = spacing * 0.25
        ax.plot([x1, x2], [y, y], color='black', linewidth=0.6)
        ax.plot([x1, x1], [y - tick_height, y], color='black', linewidth=0.6)
        ax.plot([x2, x2], [y - tick_height, y], color='black', linewidth=0.6)
        ax.text((x1 + x2) / 2, y + tick_height * 0.2, text_label,
                ha='center', va='bottom', fontsize=6)

        max_y = max(max_y, y + tick_height)
        next_y += spacing

    if max_y > upper:
        ax.set_ylim(lower, max_y + spacing * 0.2)