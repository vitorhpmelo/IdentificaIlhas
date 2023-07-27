import numpy as np

def insert_column(matrix, new_column, position):
    num_rows, num_cols = matrix.shape

    if len(new_column) != num_rows:
        raise ValueError("New column length must match the number of rows in the matrix.")

    if position < 0 or position > num_cols:
        raise ValueError("Invalid position. Position should be between 0 and the number of columns in the matrix.")

    left_half = matrix[:, :position]
    right_half = matrix[:, position:]

    new_matrix = np.concatenate((left_half, new_column[:, np.newaxis], right_half), axis=1)

    return new_matrix

