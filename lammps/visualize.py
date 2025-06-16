"""
Visualize LAMMPS trajectory from dump.springs.lammpstrj
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def read_lammps_dump(filename):
    """Read LAMMPS trajectory file"""
    frames = []
    with open(filename) as f:
        while True:
            # Read frame header
            header = f.readline()
            if not header:
                break
                
            # Skip next 8 lines (metadata)
            for _ in range(8):
                f.readline()
                
            # Read atom coordinates
            coords = []
            while True:
                line = f.readline().strip()
                if not line or line.startswith('ITEM:'):
                    break
                parts = line.split()
                coords.append([float(parts[2]), float(parts[3]), float(parts[4])])
                
            frames.append(np.array(coords))
    return frames

def visualize_frames(frames):
    """Plot 3D visualization of trajectory"""
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot first frame
    frame = frames[0]
    ax.scatter(frame[:,0], frame[:,1], frame[:,2], c='b', marker='o')
    
    # Plot last frame
    frame = frames[-1]
    ax.scatter(frame[:,0], frame[:,1], frame[:,2], c='r', marker='o')
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('LAMMPS Trajectory (Blue=Initial, Red=Final)')
    plt.savefig('trajectory_plot.png')
    print("Saved visualization to trajectory_plot.png")

if __name__ == "__main__":
    frames = read_lammps_dump("dump.springs.lammpstrj")
    print(f"Read {len(frames)} frames")
    visualize_frames(frames)
