A = [
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9]
    ] # 3x3
B = [
        [10, 11, 12],
        [13, 14, 15],
        [16, 17, 18]
    ] # 3x3

# helper function to print matrices
def print_matrix(matrix):
    for i, row in enumerate(matrix):
        print(row, "" if i < len(matrix) - 1 else "\n")

# print test matrices
print("Matrix A:")
print_matrix(A)

print("Matrix B:")
print_matrix(B)


# define functions

def dot_product(matrix):
    result = [None] * len(A[0])
    for i, vector in enumerate(matrix):
        for j, value_at_ij in enumerate(vector):
            result[j] = value_at_ij if i == 0 else result[j] * value_at_ij
    return sum(result)

def cross_product(A, B):
    Ax, Ay, Az = A
    Bx, By, Bz = B
    return [Ay*Bz - Az*By, Az*Bx - Ax*Bz, Ax*By - Ay*Bx]

def multiply_matrix(A, B):
    # manual... no loops to understand what is going where
    # make things a bit easier to see...
    # first matrix A
    # a1, a2, a3 = A[0]
    # b1, b2, b3 = A[1]
    # c1, c2, c3 = A[2]
    #
    # second matrix B
    # d1, d2, d3 = B[0]
    # e1, e2, e3 = B[1]
    # f1, f2, f3 = B[2]
    #
    # result = [
    #     [(a1*d1 + a2*e1 + a3*f1), (a1*d2 + a2*e2 + a3*f2), (a1*d3 + a2*e3 + a3*f3)],
    #     [(b1*d1 + b2*e1 + b3*f1), (b1*d2 + b2*e2 + b3*f2), (b1*d3 + b2*e3 + b3*f3)],
    #     [(c1*d1 + c2*e1 + c3*f1), (c1*d2 + c2*e2 + c3*f2), (c1*d3 + c2*e3 + c3*f3)]
    # ]

    # Get the dimensions of the matrices
    m, n, k, l = len(A), len(A[0]), len(B), len(B[0])

    # Check if the matrices can be multiplied
    if k != n:
        raise ValueError("Number of columns in A must be equal to number of rows in B")

    # Initialize the result matrix with None
    result = [[None for _ in range(l)] for _ in range(m)]

    # Perform matrix multiplication using loops for any practical size of matrices
    for i in range(m):
        for j in range(l):
            for o in range(k):
                a, b = A[i][o], B[o][j]
                # if either values of A r B are not numbers, raise an error
                if not (isinstance(a, (int, float)) and isinstance(b, (int, float))):
                    raise ValueError("All values in A and B must be numbers")
                ab = a * b
                result[i][j] = ab if result[i][j] is None else result[i][j] + ab

    return result

# print results

print("dot product of A:")
print(dot_product(A), "\n")

print("dot product of B:")
print(dot_product(B), "\n")

print("cross product of [1, 2, 3] and [4, 5, 6]:")
print(cross_product([1, 2, 3], [4, 5, 6]), "\n")

print("multiply matrix A and B:")
print_matrix(multiply_matrix(A, B))
