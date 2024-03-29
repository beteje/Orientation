#!/usr/bin/env python3
import numpy as np
from vispy import app, gloo
from vispy.util.transforms import perspective, translate, rotate
from matplotlib import pyplot as plt
from matplotlib.colors import hsv_to_rgb
import imageio
from pygifsicle import optimize

vertex = """
// Uniforms
// ------------------------------------
uniform mat4 u_model;
uniform mat4 u_view;
uniform mat4 u_projection;
uniform float u_size;

// Attributes
// ------------------------------------
attribute vec3  a_position;
attribute vec4  a_color;
attribute float a_size;

// Varyings
// ------------------------------------
varying vec4  v_color;
varying float v_size;


// Main 
// ------------------------------------
void main (void)
{
    v_size = a_size * u_size;
    v_color = a_color;
    gl_Position = u_projection * u_view * u_model * vec4(a_position,1.0);
    gl_PointSize = v_size;
}
"""

fragment = """
// Varyings 
// ------------------------------------
varying vec4 v_color;
varying float v_size;

// Main 
// ------------------------------------
void main()
{
    if (v_size >= 10000) {
        gl_FragColor = v_color;
    } else {
        float r = length((gl_PointCoord.xy - vec2(0.5, 0.5))*v_size);
        r -= v_size/2;
        if ( r > 0) {
            discard;
        } else {
            gl_FragColor = v_color;
        }
    }
} """
            
# -----------------------------------------------------------------------------
def axes():
    # define the vertices
    V = np.zeros(6, [('a_position', np.float32, 3),
              ('a_color', np.float32, 4),
              ('a_size', np.float32)])
    # set the vertex positions (origin is repeated with each colour)    
    V['a_position'] = [[-1, -1, -1], [1, -1, -1], 
                       [-1, -1, -1], [-1, 1, -1], 
                       [-1, -1, -1], [-1, -1, 1]]
    # set the vertex colours red for x-axis, green for y-axis and blue for z-axis
    V['a_color'] = [[1, 0, 0, 1], [1, 0, 0, 1], 
                    [0, 0.7, 0, 1], [0, 0.7, 0, 1], 
                    [0, 0, 1, 1], [0, 0, 1, 1]]
    V['a_size'] = np.ones(6) * 100000
    
    # set the vertex indices for the connecting lines
    I = np.arange(6, dtype=np.uint32)

    return V, I

# -----------------------------------------------------------------------------
def plot_color_wheel(az_lim, ev_lim):
    # define the azimuth and elevation values
    az_len = az_lim[1]-az_lim[0]
    ev_len = ev_lim[1]-az_lim[0]
    azth= np.linspace(0., 1., az_len)
    elvtn = np.linspace(0., 1., ev_len)
    # convert the angles to rgb colors
    rgb_array = np.zeros((az_len, ev_len, 3))
    for i in range(az_len):
        for j in range(ev_len):
            if elvtn[j] <= 0.5:
                rgb_array[i, j, :] = hsv_to_rgb([azth[i], 2*elvtn[j], 1.0])
            else:
                rgb_array[i, j, :] = hsv_to_rgb([azth[i], 1.0, 1.5 - elvtn[j]])
    
    # create figure and set axes properties
    fig, ax = plt.subplots()
    ax.set_xlim((0, ev_len))
    ax.set_xticks(np.linspace(0, ev_len, num=4))
    ax.set_xticklabels(np.linspace(ev_lim[0], ev_lim[1], num=4).astype(np.int32))
    
    ax.set_ylim((0, az_len))
    ax.set_yticks(np.linspace(0, az_len, num=4))
    ax.set_yticklabels(np.linspace(az_lim[0], az_lim[1], num=4).astype(np.int32))
    
    ax.set_xlabel('Elevation', fontsize=30, labelpad=0, color='k')
    ax.set_ylabel('Azimuth', fontsize=30, labelpad=0, color='k')
    
    ax.tick_params(direction='out', length=2, width=0, labelsize=20, pad=0, colors='k')
    
    ax.imshow(rgb_array)
    plt.show()

# -----------------------------------------------------------------------------
class Canvas(app.Canvas):
    def __init__(self):
        app.Canvas.__init__(self, keys='interactive', size=(600, 600))
        
        # set the vertices and indices for the axes lines
        self.vAxes, self.iAxes = axes()
        self.axes_buf = gloo.IndexBuffer(self.iAxes)
        
        # create 2 different programs one for the axes one for the data
        self.data_prog = gloo.Program(vertex, fragment)
        self.axes_prog = gloo.Program(vertex, fragment)
        self.axes_prog.bind(gloo.VertexBuffer(self.vAxes))
        
        # initialize the viewing position & angles
        self.translate = 5
        self.theta = -90
        self.phi = 0
        self.viewx = 0
        self.viewy = 0
        
        # set the initial view, model & projection
        self.view = translate((0, 0, -self.translate))
        self.model = np.dot(rotate(self.theta, (1, 0, 0)),
                           rotate(self.phi, (0, 1, 0)))
        self.projection = np.eye(4, dtype=np.float32)
        self.projection = np.eye(4, dtype=np.float32)
        self.rescale()
        
        # add the view, model & size to each program
        self.data_prog['u_view'] = self.view
        self.axes_prog['u_view'] = self.view
        self.data_prog['u_model'] = self.model
        self.axes_prog['u_model'] = self.model
        self.data_prog['u_size'] = 5 / self.translate
        self.axes_prog['u_size'] = 5 / self.translate
        
        gloo.set_state('translucent', clear_color='white')
        self.show()
    
    # ---------------------------------
    def on_draw(self, event):
        gloo.clear()
        # draw the data points
        self.data_prog.draw('points')
        
        # draw the axes lines
        gloo.set_state(blend=True, depth_test=True, polygon_offset_fill=False)
        gloo.set_depth_mask(False)
        self.axes_prog.draw('lines', self.axes_buf)
        gloo.set_depth_mask(True)
    
    # ---------------------------------
    def on_resize(self, event):
        self.rescale()
        
    # ---------------------------------
    def rescale(self):
        # adjust the viewport and projection matrix based on current physical size
        gloo.set_viewport(0, 0, self.physical_size[0], self.physical_size[1])
        self.projection = perspective(45.0, self.size[0] /
                                      float(self.size[1]), 1.0, 1000.0)
        self.data_prog['u_projection'] = self.projection
        self.axes_prog['u_projection'] = self.projection
     
    # ---------------------------------
    def on_key_press(self, event):
        # rotate both the data and axes based on direction keys
        if event.key == 'Left':
            self.phi += .5
        elif event.key == 'Right':
            self.phi -= .5
        elif event.key =='Up':
            self.theta += .5
        elif event.key == 'Down':
            self.theta -= .5
        
        # adjust the model based on a rotation of theta round the x-axis and phi round the z-axis
        self.model = np.dot(rotate(self.theta, (1, 0, 0)),
                           rotate(self.phi, (0, 1, 0)))
        
        self.data_prog['u_model'] = self.model
        self.axes_prog['u_model'] = self.model
        self.update()
        
    # ---------------------------------
    def on_mouse_wheel(self, event):
        # use the mouse wheel to adjust the zoom on the view
        self.translate -= event.delta[1]/10
        self.translate = max(0.01, self.translate)
        self.view = translate((self.viewx, self.viewy, -self.translate))
        
        self.data_prog['u_view'] = self.view
        self.axes_prog['u_view'] = self.view
        # adjust the sizes of the data points
        self.data_prog['u_size'] = 5 / self.translate
        self.update()
        
    # ---------------------------------
    def on_mouse_move(self, event):
        if event.is_dragging:
            # adjust the x and y positions (of the view)
            x0, y0 = event.press_event.pos
            x1, y1 = event.last_event.pos
            dx, dy = x1 - x0, y1 - y0
            
            self.viewx += dx / 10000
            self.viewy -= dy / 10000

            self.view = translate((self.viewx, self.viewy, -self.translate))
            self.data_prog['u_view'] = self.view
            self.axes_prog['u_view'] = self.view
            self.update()
      
    # ---------------------------------
    def create_animation(self, fname):
        # sleep briefly to prevent non-exposed window warnings
        self.app.sleep(0.05)
        # Set the frames per second
        fps = 20
        # Set the number of steps and number of degrees to update by
        # (Gives approx 1 rotation every 10 seconds)
        n_steps = fps*10
        n_degrees = 360/n_steps
        writer = imageio.get_writer(fname, fps=fps)
        for i in range(n_steps * 2):
            im = self.render()      
            im = im[50:550,50:550,0:3]
            writer.append_data(im)
            if i >= n_steps:
                self.phi += n_degrees
            else:
                self.theta += n_degrees
        
            # adjust the model based on a rotation of theta round the x-axis and phi round the z-axis
            self.model = np.dot(rotate(self.theta, (1, 0, 0)),
                               rotate(self.phi, (0, 1, 0)))
            
            self.data_prog['u_model'] = self.model
            self.axes_prog['u_model'] = self.model
            self.update()
            
        writer.close() 
        optimize(fname)
      
# -----------------------------------------------------------------------------            
if __name__ == '__main__':
    # create a canvas
    c = Canvas()
    # produce some data to put into the canvas
    n = 500
    ps = c.pixel_scale
    data = np.zeros(n, [('a_position', np.float32, 3),
                        ('a_color', np.float32, 4),
                        ('a_size', np.float32)])
    tmp = 0.45 * np.random.randn(n,3)
    data['a_position'] = 2 * ((tmp - np.amin(tmp)) / (np.amax(tmp) - np.amin(tmp))) - 1
    data['a_color'] = np.concatenate((np.random.uniform(0.5, 1.00, (n, 3)), np.ones((n,1))), axis=1)
    data['a_size']   = np.ones(n) * 20. * ps
    c.data_prog.bind(gloo.VertexBuffer(data))
    c.app.run()
    
    # create a gif of the canvas rotating
    c.create_animation('bubbles.gif')