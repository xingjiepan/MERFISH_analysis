#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, CheckButtons


class ImagePlotter:
    '''This class visualize an image with multiple color channels and
    and allows users to navigate through different sections.
    '''

    def __init__(self, z_scale=10):
        self.channels = {}
        self.z = 0
        self.y = 0
        self.z_scale = z_scale
        
    def init_plot(self, figsize=None):
        stack_size = self.get_image_stack_size()
        self.fig = plt.figure(figsize=figsize, constrained_layout=True)
        widths = [1, 0.05]
        heights = [1, 1, 0.2]
        gs = self.fig.add_gridspec(nrows=3, ncols=2, width_ratios=widths, height_ratios=heights) 

        # Initialize plotting areas
        self.ax1 = self.fig.add_subplot(gs[0, 0])
        self.ax2 = self.fig.add_subplot(gs[1, 0], sharex=self.ax1)
        self.ax3 = self.fig.add_subplot(gs[0, 1])
        self.ax4 = self.fig.add_subplot(gs[1, 1])
        self.ax5 = self.fig.add_subplot(gs[2, :])

        self.ax1.set_aspect('equal')
        self.ax2.set_aspect(self.z_scale)

        # Initialize sliders
        self.slider1 = Slider(self.ax3, 'Z', 0, stack_size[0], valinit=0, orientation='vertical')
        self.slider2 = Slider(self.ax4, 'Y', 0, stack_size[2], valinit=0, orientation='vertical')

        self.slider1.on_changed(self.update_slider1)
        self.slider2.on_changed(self.update_slider2)

        # Initialize check buttons
        visibility = [self.channels[k][2] for k in self.channels]
        self.checkbox = CheckButtons(self.ax5, self.channels.keys(), visibility) 
        self.checkbox.on_clicked(self.update_checkbox)

    def update_slider1(self, val):
        self.z = int(self.slider1.val)
        self.plot_x_y_plane(preserve_image_size=True)

    def update_slider2(self, val):
        self.y = int(self.slider2.val)
        self.plot_x_z_plane(preserve_image_size=True)

    def update_checkbox(self, label):
        self.channels[label][2] = not self.channels[label][2]
        self.plot_x_y_plane(preserve_image_size=True)
        self.plot_x_z_plane(preserve_image_size=True)

    def add_image_stack(self, channel_name, image_stack, color):
        # Normalize the RGB color
        if not (type(color) is str or color is None):
            color = np.array(color) / max(color) 
        
        self.channels[channel_name] = [image_stack, color, True] 

    def get_image_stack_size(self):
        '''Get the size of the image stack.
        Return a tuple of (n_z, n_x, n_y, n_channels).
        '''
        n_channels = len(self.channels.keys())

        Zs = []
        Xs = []
        Ys = []
        for channel_name in self.channels:
            image_stack, color_map, visible = self.channels[channel_name]
            Zs.append(image_stack.shape[0])
            Xs.append(image_stack.shape[1])
            Ys.append(image_stack.shape[2])

        if n_channels == 0:
            return (0, 0, 0, 0)

        if max(Zs) != min(Zs) or max(Xs) != min(Xs) or max(Ys) != min(Ys):
            raise Exception('The image sizes of different channels do not match.')

        return (Zs[0], Xs[0], Ys[0], n_channels)

    def n_active_channels(self):
        n_active = 0
        for channel_name in self.channels:
            n_active += self.channels[channel_name][2]
        return n_active

    def get_colors_for_clip(self, clip, color):
        # Use the defined channel color
        if type(color) is np.ndarray: 
            return clip[:, :, np.newaxis] * color[np.newaxis, np.newaxis, :]
       
        # Use the chaotic color generator
        if color == 'chaotic':
            def chaotic_func(x, shift):
                return np.sin(100 / (x + 1e-20 + shift * 0.01)) * (x > 0)

            return np.stack((chaotic_func(clip, 0), chaotic_func(clip, 1), chaotic_func(clip, 2)), axis=-1)

        # Default color
        return np.stack((clip, clip, clip), axis=-1)

    def plot_x_y_plane(self, preserve_image_size=False):
        # Get the current image size
        if preserve_image_size:
            x_lim = self.ax1.get_xlim()
            y_lim = self.ax1.get_ylim()

        # Update the plot 
        self.ax1.cla()
        stack_size = self.get_image_stack_size()
        img = np.zeros((stack_size[2], stack_size[1], 3))
        
        for channel_name in self.channels:
            image_stack, color, visible = self.channels[channel_name]
            
            if visible:
                v_max = np.max(image_stack)
                v_min = np.min(image_stack)
                clip = np.transpose(image_stack[self.z] - v_min) / (v_max - v_min + 1e-20)
                img += self.get_colors_for_clip(clip, color)
        
        self.ax1.imshow(img / (self.n_active_channels() + 1e-20))

        # Reset the image size
        if preserve_image_size:
            self.ax1.set_xlim(*x_lim)
            self.ax1.set_ylim(*y_lim)

        self.ax1.set_xlabel('X')
        self.ax1.set_ylabel('Y')

    def plot_x_z_plane(self, preserve_image_size=False):
        # Get the current image size
        if preserve_image_size:
            x_lim = self.ax2.get_xlim()
            y_lim = self.ax2.get_ylim()

        # Update the plot 
        self.ax2.cla()
        stack_size = self.get_image_stack_size()
        
        img = np.zeros((stack_size[0], stack_size[1], 3))
        
        for channel_name in self.channels:
            image_stack, color, visible = self.channels[channel_name]
            
            if visible:
                v_max = np.max(image_stack)
                v_min = np.min(image_stack)
                clip = (image_stack[:, :, self.y] - v_min) / (v_max - v_min + 1e-20)
                img += self.get_colors_for_clip(clip, color)
        
        self.ax2.imshow(img / (self.n_active_channels() + 1e-20))

        # Reset the image size
        if preserve_image_size:
            self.ax2.set_xlim(*x_lim)
            self.ax2.set_ylim(*y_lim)
        
        self.ax2.set_xlabel('X')
        self.ax2.set_ylabel('Z')
        self.ax2.set_aspect(self.z_scale)

    def show(self, figsize=None):
        self.init_plot(figsize=figsize)
        self.plot_x_y_plane()
        self.plot_x_z_plane()


