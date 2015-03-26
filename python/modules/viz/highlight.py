import sys
import matplotlib.patches as patches

def highlight_x(ax, highlight_range, highlight_color='magenta', label=None):
    """Highlights an x-range on a given axes object"""

    rect = patches.Rectangle((highlight_range[0], 0), highlight_range[1] - highlight_range[0], ax.get_ylim()[1], facecolor=highlight_color, edgecolor='none', alpha=0.5)
    ax.add_patch(rect)

    if label is not None:
        ax.text(highlight_range[0] + ((highlight_range[1] - highlight_range[0]) / 10), 0.9 * ax.get_ylim()[1], label, rotation=90, color=highlight_color)

def highlight_y(ax, highlight_range, highlight_color='magenta', label=None):
    """Highlights an y-range on a given axes object"""

    rect = patches.Rectangle((0, highlight_range[0]), ax.get_xlim()[1], highlight_range[1] - highlight_range[0], facecolor=highlight_color, edgecolor='none', alpha=0.5)
    ax.add_patch(rect)

    if label is not None:
        ax.text(0.9 * ax.get_xlim()[1], highlight_range[0] + ((highlight_range[1] - highlight_range[0]) / 2), label, color=highlight_color)

