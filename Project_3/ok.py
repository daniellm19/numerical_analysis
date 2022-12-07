# Prompt the user for the size of the matrix
size = int(input("Enter the size of the matrix: "))

# Create an empty matrix with the specified size
matrix = []
for i in range(size):
    matrix.append([0] * size)

# Fill in the matrix with the values described in the question
for i in range(size):
    for j in range(size):
        if i == 0:
            matrix[i][j] = 4
        elif i == size - 1:
            matrix[i][j] = 2
        elif j == 0:
            matrix[i][j] = 1
        elif j == size - 1:
            matrix[i][j] = 3
        else:
            matrix[i][j] = 5

# Print the matrix to the console
for row in matrix:
    print(row)

