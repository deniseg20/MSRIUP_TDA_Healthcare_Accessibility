import numpy as np
import pandas as pd

def extend_array_with_nan_diagonal(output_file=None):
    """
    Read a 301x301 npy file and extend it to 411x411 with NaN diagonal
    
    Parameters:
    input_file (str): Path to input .npy file
    output_file (str, optional): Path to save extended array
    
    Returns:
    numpy.ndarray: Extended 411x411 array
    """
    
    # Read the original 301x301 array
    original_array = np.load('persistent_homology/distance_matrices/fqhc_d_matrix.npy')
    
    # Verify dimensions
    if original_array.shape != (301, 301):
        raise ValueError(f"Expected 301x301 array, got {original_array.shape}")
    
    # Create new 411x411 array (301 + 110 = 411) filled with large number
    extended_array = np.full((411, 411), 10000.0)
    
    # Copy original data to top-left corner (indices 0:301, 0:301)
    extended_array[:301, :301] = original_array
    
    # Set the entire 110x110 square (indices 301:411, 301:411) to NaN
    # This prevents all connections within the new region
    extended_array[301:411, 301:411] = np.nan
    
    # Set the diagonal of the new region to 0 (distance from each point to itself)
    for i in range(110):
        extended_array[301 + i, 301 + i] = 0.0
    
    # Optionally save the result
    if output_file:
        np.save(output_file, extended_array)
        print(f"Extended array saved to {output_file}")
    
    return extended_array

# Example usage
if __name__ == "__main__":
    # Load and extend the array
    input_filename = "fqhc_d_matrix.npy"  # Replace with your file path
    output_filename = "extended_array.npy"  # Optional output file
    
    try:
        result = extend_array_with_nan_diagonal(output_filename)
        
        print(f"Original shape: (301, 301)")
        print(f"Extended shape: {result.shape}")
        print(f"Values at original region (0:301, 0:301): preserved")
        print(f"Values in extended square (301:411, 301:411): NaN except diagonal")
        print(f"Diagonal values in extended region: 0")
        print(f"Values in other extended regions: 10000")
        
        # Show some sample values
        print(f"\nSample original value at [0,0]: {result[0,0]}")
        print(f"Sample extended value at [301,250]: {result[301,250]} (should be 10000)")
        print(f"Sample extended value at [301,302]: {result[301,302]} (should be NaN)")
        print(f"Sample diagonal value at [301,301]: {result[301,301]} (should be 0)")
        print(f"Sample diagonal value at [350,350]: {result[350,350]} (should be 0)")
        
    except FileNotFoundError:
        print(f"File {input_filename} not found. Please update the filename.")
    except ValueError as e:
        print(f"Error: {e}")
    