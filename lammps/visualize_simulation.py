"""
Visualize LAMMPS trajectory from dump.springs.lammpstrj using seaborn (2D)
"""
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import glob

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
                coords.append([float(parts[2]), float(parts[3])])  # Only x,y for 2D
                
            frames.append(np.array(coords))
    return frames

def visualize_frames(frames, step_num):
    """Plot 2D visualization of trajectory using seaborn"""
    plt.figure(figsize=(10, 8))
    
    # Plot points
    frame = frames[-1]
    sns.scatterplot(x=frame[:,0], y=frame[:,1], color='red', s=100)
    
    # Simple line connections (for visualization purposes)
    plt.plot(frame[:,0], frame[:,1], 'b-', linewidth=1)
    
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(f'LAMMPS Simulation - Step {step_num}')
    plt.savefig(f"step_{step_num}.png")
    plt.close()

if __name__ == "__main__":
    # Process all dump files from the loop
    dump_files = sorted(glob.glob("dump.springs_loop.*.lammpstrj"))
    
    for i, dump_file in enumerate(dump_files):
        frames = read_lammps_dump(dump_file)
        visualize_frames(frames, i)
        print(f"Visualized step {i} - saved to step_{i}.png")
