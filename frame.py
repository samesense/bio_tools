#!/usr/bin/env python
"""
Frame classes for customizing frame borders that surround the plot axes.
"""
import numpy as np

import matplotlib.axes as maxes
import matplotlib.pyplot as plt
import matplotlib.artist as martist
import matplotlib.collections as col
import matplotlib.projections as projections


class Frame(martist.Artist):
    """Draw frame along the edges of the axes patch.
    
    Frame position can be controlled upon initialization or by setting
    `positions` property with a list of positions
        ['left', 'right', 'top', 'bottom' | 'all']
    """
    _position_list = ('left', 'right', 'top', 'bottom')
    def __init__(self, axes, positions=('left', 'bottom'), **kwargs):
        """
        `positions` is a list of strings of frame positions to plot.
            ['left', 'right', 'top', 'bottom' | 'all']
        """
        super(Frame, self).__init__()
        # TODO: allow more keyword configuration
        self.axes = axes
        
        rc = plt.rcParams
        self.color = kwargs.pop('color', rc['axes.edgecolor'])
        self.linewidth = kwargs.pop('linewidth', rc['axes.linewidth'])
        self.linestyle = kwargs.pop('linestyle', 'solid')
        self.positions = positions

    def get_data(self):
        """Convenience method returns tuple of (x, y) data in `self.axes`"""
        x, y = [], []
        ax = self.axes
        for artist in (ax.lines, ax.patches):
            if not artist == []:
                x.append(np.concatenate([a.get_xdata() for a in artist]))
                y.append(np.concatenate([a.get_ydata() for a in artist]))
        # TODO: get scatter data from ax.collections
        return (np.concatenate(x), np.concatenate(y))

    def _set_frame_position(self, positions):
        """Set positions where frame will be drawn.

        `positions` is a list of strings of frame positions to plot.
            ['left', 'right', 'top', 'bottom' | 'all']
        """
        self._frame_on = self._frame_dict_from(positions)

    def _get_frame_position(self):
        return [p for p in self._position_list if self._frame_on[p]]

    # xposition tuples turn on frame for (bottom, top)
    _xposition_pairs = {(True, False): 'bottom', (False, True): 'top',
                        (True, True): 'both', (False, False): 'none'}
    def _get_xposition(self, frame_on=None):
        """Returns position that matches `XAxis.set_ticks_position` inputs.

        `frame_on` is a dict that matches frame positions with bools.
        """
        if frame_on is None:
            frame_on = self._frame_on
        return self._xposition_pairs[(frame_on['bottom'], frame_on['top'])]

    # yposition tuples turn on frame for (left, right)
    _yposition_pairs = {(True, False): 'left', (False, True): 'right',
                        (True, True): 'both', (False, False): 'none'}
    def _get_yposition(self, frame_on=None):
        """Returns position that matches `YAxis.set_ticks_position` inputs.

        `frame_on` is a dict that matches frame positions with bools.
        """
        if frame_on is None:
            frame_on = self._frame_on
        return self._yposition_pairs[(frame_on['left'], frame_on['right'])]

    def _frame_dict_from(self, positions):
        """Parse `positions` and return xposition, yposition tuple

        `positions` is a list of strings of frame positions to plot.
            ['left', 'right', 'top', 'bottom' | 'all']
        """
        frame_dict = dict.fromkeys(self._position_list, False)
        
        if 'all' in positions:
            frame_dict = dict.fromkeys(self._position_list, True)
        else:
            for position in positions:
                frame_dict[position] = True
        return frame_dict

    def _set_ticks(self):
        """Overide this method to customize tick positioning."""
        # Draw ticks on axes only where a frame is drawn
        self.axes.xaxis.set_ticks_position(self._get_xposition())
        self.axes.yaxis.set_ticks_position(self._get_yposition())

    _frame_lines = dict(bottom=[(0., 0.), (1., 0.)], top=[(0., 1.), (1., 1.)],
                        left=[(0., 0.), (0., 1.)], right=[(1., 0.), (1., 1.)])
    def _make_frame(self):
        """Get axis frame specified by `self._frame_on`."""
        lines = [self._frame_lines[p] for p in self._position_list
                                      if self._frame_on[p]]
        frame_lines = col.LineCollection(segments=lines,
                                         linewidths=[self.linewidth],
                                         colors=[self.color])
        frame_lines.set_transform(self.axes.transAxes)
        return frame_lines

    def draw(self, renderer):
        if not self.get_visible():
            return
        self._set_ticks()
        frame = self._make_frame()
        frame.draw(renderer)

    positions = property(_get_frame_position, _set_frame_position)


class FrameAxes(maxes.Axes):
    """Override Axes so that `frame` attribute is a custom FrameArtist
    """
    name = 'frameaxes'

    def cla(self):
        """Override default Axes.frame and clear patch edgecolor"""
        if not hasattr(self, '_frame'):
            self._frame = Frame(self)
        super(FrameAxes, self).cla()
        self.patch.set_edgecolor('none')
        # retrieve `_frame` because Axes.cla() overwrites `frame`
        self.frame = self._frame
        # self.frame.set_zorder(10)

    def set_frame(self, frame_artist, **kwargs):
        """Set artist to draw axes frame.

        `frame_artist` should be subclass of `FrameArtist`.
        """
        # save frame artist as `_frame` because Axes.cla() overwrites `frame`
        self.frame = self._frame = frame_artist(self, **kwargs)

    def get_frame(self):
        return self.frame
projections.register_projection(FrameAxes)


if __name__ == '__main__':
    # a quick example on different ways to draw the frame
    from numpy.random import rand
    
    positions_test = (None, ['all'], ['top', 'right'])
    names = ('L frame (default)', 'normal frame', 'top, right')
    sub_pos = (221, 223, 122)

    x = rand(10)
    y = rand(10)
    for positions, name, sub in zip(positions_test, names, sub_pos):
        ax = plt.subplot(sub, projection='frameaxes')
        if positions is not None:
            ax.frame.positions = positions
        ax.plot(x, y, 'o')
        ax.set_ylabel(name)
    fig = plt.gcf()
    # remove background color to make the frames stand out
    fig.set_facecolor('w')
    fig.subplots_adjust(hspace=0.5, wspace=0.5)
    plt.show()
