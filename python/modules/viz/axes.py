def set_spines(ax, visible=False, spine_set=None):
    
    if spine_set is None:
        spine_set = ['left', 'right', 'top', 'bottom']
    for s in spine_set:
        ax.spines[s].set_visible(visible)

def no_ticks(ax, no_x=True, no_y=True):
    if no_x:
        ax.xaxis.set_ticks_position('none')
    if no_y:
        ax.yaxis.set_ticks_position('none')

def clean_axis(ax, right=True, left=False, top=True, bottom=False, allx=False):
    if right or allx:
        ax.spines['right'].set_visible(False)
    if left or allx:
        ax.spines['left'].set_visible(False)
        ax.set_yticks([])
    if top or allx:
        ax.spines['top'].set_visible(False)
    if bottom or allx:
        ax.spines['bottom'].set_visible(False)
        ax.set_xticks([])
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

def set_ticks_outer(ax, x=True, y=True):
    if y:
        ax.get_yaxis().set_tick_params(which='both', direction='out')
    if x:
        ax.get_xaxis().set_tick_params(which='both', direction='out')

def set_common_range(axes, x=False, y=True, setmax=True, setmin=True):

    if y:
        y_max = None 
        y_min = None
        for ax in axes:
            y_min = ax.get_ylim()[0] if y_min is None else min(y_min, ax.get_ylim()[0])
            y_max = max(y_max, ax.get_ylim()[1])
        for ax in axes:
            if setmin and setmax:
                ax.set_ylim([y_min, y_max])
            elif setmin:
                ax.set_ylim([y_min, ax.get_ylim()[1]])
            elif setmax:
                ax.set_ylim([ax.get_ylim()[0], y_max])
    if x:        
        x_max = None 
        x_min = None
        for ax in axes:
            x_min = ax.get_xlim()[0] if x_min is None else min(x_min, ax.get_xlim()[0])
            x_max = max(x_max, ax.get_xlim()[1])
        for ax in axes:
            if setmin and setmax:
                ax.set_xlim([x_min, x_max])
            elif setmin:
                ax.set_xlim([x_min, ax.get_xlim()[1]])
            elif setmax:
                ax.set_xlim([ax.get_xlim()[0], x_max])
         
def label_bars(ax, rects, labels, rotation=0, fontsize=10):
    
    assert len(labels) == len(rects)

    for i, r in enumerate(rects):
        #height = r.get_height()
        height = ax.get_ylim()[1] * 0.5
        ax.text(r.get_x() + r.get_width() / 2, 
                1.05 * height,
                '%s' % str(labels[i]),
                ha='center',
                va='bottom',
                rotation=rotation,
                fontsize=fontsize )

