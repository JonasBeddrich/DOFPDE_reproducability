import pyvista as pv
import numpy as np 
import glob 
from pyvista import examples
from tqdm import tqdm
import time
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd 

def colormap(): 
    cooltowarmcolors = [
        (0.000000, 0.000000, 0.349020), 
        (0.039216, 0.062745, 0.380392), 
        (0.062745, 0.117647, 0.411765), 
        (0.090196, 0.184314, 0.450980), 
        (0.125490, 0.262745, 0.501961), 
        (0.160784, 0.337255, 0.541176), 
        (0.200000, 0.396078, 0.568627),      
        (0.239216, 0.454902, 0.600000),
        (0.286275, 0.521569, 0.650980),
        (0.337255, 0.592157, 0.701961),
        (0.388235, 0.654902, 0.749020),
        (0.466667, 0.737255, 0.819608),
        (0.572549, 0.819608, 0.878431),
        (0.654902, 0.866667, 0.909804),
        (0.752941, 0.917647, 0.941176),
        (0.823529, 0.956863, 0.968627),
        (0.941176, 0.984314, 0.988235),
        (0.988235, 0.960784, 0.901961),
        (0.988235, 0.945098, 0.850980),
        (0.980392, 0.898039, 0.784314),
        (0.968627, 0.835294, 0.698039),
        (0.949020, 0.733333, 0.588235),
        (0.929412, 0.650980, 0.509804),
        (0.909804, 0.564706, 0.435294),
        (0.878431, 0.458824, 0.352941),
        (0.839216, 0.388235, 0.286275),
        (0.760784, 0.294118, 0.211765),
        (0.701961, 0.211765, 0.168627),
        (0.650980, 0.156863, 0.129412),
        (0.600000, 0.094118, 0.094118),
        (0.549020, 0.066667, 0.098039),
        (0.501961, 0.050980, 0.125490),
        (0.450980, 0.052902, 0.172549),
        (0.400000, 0.054902, 0.192157),
        (0.349020, 0.070588, 0.211765)
    ]

    cooltowarmcmap = LinearSegmentedColormap.from_list("cooltowarm", cooltowarmcolors)

    return cooltowarmcmap

def scaled_colormap(): 
    cmap = colormap() 
    x = np.linspace(-5,5,256)
    sigmoid = 1 / (1 + np.exp(-x))
    # normalize 
    sigmoid -= sigmoid[0]
    sigmoid /= sigmoid[-1]
    # get colors 
    colors = cmap(sigmoid)
    scaled_colormap = LinearSegmentedColormap.from_list("scaledcooltowarm", colors)
    return scaled_colormap

def main():

    # Create sample colormaps for demo
    cooltowarmcmap = scaled_colormap() 
    
    filelist = glob.glob('data/*1000*/*.pvd')
    filelist.sort()
    print(filelist)

    for filepath in filelist: 
        start = time.time()    
        filename = filepath.split("/")[-1]
        reader = pv.get_reader(filepath)
        print(reader.time_values)
        for t in reader.time_values:
            reader.set_active_time_value(t)

            print("reading mesh")
            print(time.time() - start)
            mesh = reader.read()

            print("scaling")
            print(time.time() - start)
            scaled_blocks = []
            for block in tqdm(mesh): 
                blockcp = block.copy()
                blockcp.scale([1,1,10])
                scaled_blocks.append(blockcp)
            
            scaled = pv.MultiBlock(scaled_blocks)
            # scaled.plot(scalars='u')

            print("surface")
            print(time.time() - start)
            combined = pv.PolyData()
            for block in tqdm(scaled): 
                combined = combined + block 

            surface = combined.extract_surface()
            # surface.plot()

            vabs = np.max([np.abs(surface["u"].max()), np.abs(surface["u"].min())])
            warpedsurface = surface.warp_by_scalar('u', factor= - 300 / vabs)

            print("plot")
            print(time.time() - start)

            scalar_bar_args = {
                "position_x": 0.2,   # move right (horizontal position)
                "position_y": 0.27,  # move up from bottom
                "width": 0.3,        # make it narrow
                "height": 0.08,        # make it short
                "color": "black", 
                "n_labels": 3,             # 👈 Exactly three tick labels
                "label_font_size": 96,
                "title_font_size": 96,
                "fmt": "%.1e"     
            }

            # plotter = pv.Plotter()
            plotter = pv.Plotter(off_screen=True, window_size=(3000, 2000))
            plotter.set_background('white', top=None)

            # plotter.add_mesh(mesh)
            plotter.add_mesh(warpedsurface, scalars='u', cmap=cooltowarmcmap, clim=(-vabs, vabs), scalar_bar_args=scalar_bar_args)

            # Set camera position: (position, focal_point, view_up)
            plotter.camera_position = [
                (18891.030680957912, -8209.117491297726, -9665.679298636129),
                (476.5681484633901, 10683.534144139054, 944.1738152586305),
                (-0.30933019099835724, 0.21828283369999812, -0.9255633081798453)
            ]

            # Render and save screenshot
            plotter.enable_anti_aliasing()
            plotter.show(auto_close=False)  # required to render before screenshot
            plotter.screenshot(filename[:-4] + "_t_" + str(t) + ".png")
            plotter.close()

            final_camera_position = plotter.camera_position
            print("Final camera position:", final_camera_position)
            print(time.time() - start)

def cut_out():
    # Create sample colormaps for demo
    cooltowarmcmap = scaled_colormap() 
    
    filelist = glob.glob('data/*/*.pvd')
    filelist.sort()
    print(filelist)

    for filepath in filelist: 
        start = time.time()    
        filename = filepath.split("/")[-1]
        reader = pv.get_reader(filepath)

        print("reading mesh")
        print(time.time() - start)
        mesh = reader.read()

        print("scaling")
        print(time.time() - start)
        scaled_blocks = []
        for block in tqdm(mesh): 
            blockcp = block.copy()
            blockcp.scale([1,1,10])
            scaled_blocks.append(blockcp)
        
        scaled = pv.MultiBlock(scaled_blocks)
        # scaled.plot(scalars='u')

        print("surface")
        print(time.time() - start)
        combined = pv.PolyData()
        for block in tqdm(scaled): 
            combined = combined + block 

        vabs = np.max([np.abs(combined["u"].max()), np.abs(combined["u"].min())])
      
        points = np.asarray([
            (3500, 2500, 0), 
            (2875, 7800, 0), 
            (2400, 11500, 0), 
            (3100, 15500, 0)
        ])

        S1 = combined.slice(normal=(0,-1,0), origin = points[0])
        S2 = combined.slice(normal=(0,-1,0), origin = points[1])
        S3 = combined.slice(normal=(0,-1,0), origin = points[2])
        S4 = combined.slice(normal=(0,-1,0), origin = points[3])

        C1 = combined.clip(normal=(1,0,0),origin=points[3]).clip(normal=(0,-1,0), origin = points[3])
        C2 = combined.clip(normal=(0,1,0),origin=points[3]).clip(normal=(0,-1,0), origin = points[2])
        C3 = combined.clip(normal=(0,1,0),origin=points[2]).clip(normal=(0,-1,0), origin = points[1])
        C4 = combined.clip(normal=(0,1,0),origin=points[1]).clip(normal=(0,-1,0), origin = points[0])
        C5 = combined.clip(normal=(1,0,0),origin=points[0]).clip(normal=(0,1,0), origin = points[0])

        tmp = points[3]-points[2]
        C2 = C2.clip(normal=(tmp[1], -tmp[0], 0), origin=points[3])

        tmp = points[2]-points[1]
        C3 = C3.clip(normal=(tmp[1], -tmp[0], 0), origin=points[2])
        
        tmp = points[1]-points[0]
        C4 = C4.clip(normal=(tmp[1], -tmp[0], 0), origin=points[1])
        
        print("plot")
        print(time.time() - start)

        scalar_bar_args = {
            "position_x": 0.17,   # move right (horizontal position)
            "position_y": 0.27,  # move up from bottom
            "width": 0.35,        # make it narrow
            "height": 0.08,        # make it short
            "color": "black"
        }

        # plotter = pv.Plotter()
        plotter = pv.Plotter(off_screen=True, window_size=(3000, 2000))
        plotter.set_background('white', top=None)

        plotter.add_mesh(C1, scalars='u', cmap=cooltowarmcmap, clim=(-vabs, vabs), scalar_bar_args=scalar_bar_args)
        plotter.add_mesh(C2, scalars='u', cmap=cooltowarmcmap, clim=(-vabs, vabs), scalar_bar_args=scalar_bar_args)
        plotter.add_mesh(C3, scalars='u', cmap=cooltowarmcmap, clim=(-vabs, vabs), scalar_bar_args=scalar_bar_args)
        plotter.add_mesh(C4, scalars='u', cmap=cooltowarmcmap, clim=(-vabs, vabs), scalar_bar_args=scalar_bar_args)
        plotter.add_mesh(C5, scalars='u', cmap=cooltowarmcmap, clim=(-vabs, vabs), scalar_bar_args=scalar_bar_args)

        plotter.add_mesh(S1, scalars='u', cmap=cooltowarmcmap, clim=(-vabs, vabs), scalar_bar_args=scalar_bar_args)
        plotter.add_mesh(S2, scalars='u', cmap=cooltowarmcmap, clim=(-vabs, vabs), scalar_bar_args=scalar_bar_args)
        plotter.add_mesh(S3, scalars='u', cmap=cooltowarmcmap, clim=(-vabs, vabs), scalar_bar_args=scalar_bar_args)
        plotter.add_mesh(S4, scalars='u', cmap=cooltowarmcmap, clim=(-vabs, vabs), scalar_bar_args=scalar_bar_args)

        plotter.camera_position = [
            (18891.030680957912, -8209.117491297726, -9665.679298636129),
            (476.5681484633901, 10683.534144139054, 944.1738152586305),
            (-0.30933019099835724, 0.21828283369999812, -0.9255633081798453)
        ]

        # Render and save screenshot
        plotter.enable_anti_aliasing()
        plotter.show(auto_close=False)  # required to render before screenshot
        plotter.screenshot(filename[:-4] + "_cut_view.png")
        plotter.close()

        # final_camera_position = plotter.camera_position
        # print("Final camera position:", final_camera_position)
        # print(time.time() - start)        

def depth_profile():
    
    start = time.time()    

    filelist = glob.glob('data/Scenario_*/*.pvd')
    print(filelist)

    for filepath in filelist[:1]: 

        filename = filepath.split("/")[-1]
        reader = pv.get_reader(filepath)

        print("reading mesh")
        print(time.time() - start)

        mesh = reader.read()

        print("scaling")
        print(time.time() - start)

        scaled_blocks = []
        for block in tqdm(mesh): 
            blockcp = block.copy()
            blockcp.scale([1,1,10])
            scaled_blocks.append(blockcp)
        scaled = pv.MultiBlock(scaled_blocks)


        flat_blocks = []
        for block in tqdm(mesh): 
            blockcp = block.copy()
            blockcp.scale([1,1,-0.01])
            flat_blocks.append(blockcp)
        flat = pv.MultiBlock(flat_blocks)

        print("surface of the domain")
        print(time.time() - start)

        combined = pv.PolyData()
        for block in tqdm(scaled): 
            combined = combined + block 
        surface = combined.extract_surface()

        surface_depth = surface.copy()
        z_coords = surface_depth.points[:, 2]
        surface_depth["depth"] = z_coords / 10

        print("water surface")
        print(time.time() - start)

        flat_combined = pv.PolyData()
        for block in tqdm(flat): 
            flat_combined = flat_combined + block 
        flat_surface = flat_combined.extract_surface()

        flat_depth = flat_surface.copy()
        z_coords = flat_depth.points[:, 2]
        flat_depth["depth"] = z_coords * -100

        print("plot")
        print(time.time() - start)

        scalar_bar_args = {
            "position_x": 0.4,   # move right (horizontal position)
            "position_y": 0.1,  # move up from bottom
            "width": 0.5,        # make it narrow
            "height": 0.1        # make it short
        }

        plotter = pv.Plotter(off_screen=True)
        plotter.add_mesh(surface_depth, scalars='depth', scalar_bar_args=scalar_bar_args)
        plotter.add_mesh(flat_depth, scalars='depth', scalar_bar_args=scalar_bar_args)

        # Set camera position: (position, focal_point, view_up)
        plotter.camera_position = [
            (11095.115571088621, 26440.891496872835, -9501.71392174851),
            (2301.2595977783203, 12050.0, 648.8344039916992),
            (-0.3410075537930717, -0.39372755402465814, -0.8536348525322903)
        ]

        # Render and save screenshot
        plotter.show(auto_close=False)  # required to render before screenshot
        # plotter.screenshot(filename[:-4] + "_depth_profile.png")
        # plotter.close()

        # print(time.time() - start)

if __name__ == "__main__":
    main()
    # cut_out()
    # depth_profile()