import numpy as np
np.set_printoptions(threshold=10000)  # Show up to 1000 elements
from numpy import matrix
from numpy import linalg

def modMatInv(A,p):       # Finds the inverse of matrix A mod p
  n=len(A)
  A=matrix(A)
  adj=np.zeros(shape=(n,n))
  for i in range(0,n):
    for j in range(0,n):
      adj[i][j]=((-1)**(i+j)*int(round(linalg.det(minor(A,j,i)))))%p
  return (modInv(int(round(linalg.det(A))),p)*adj)%p

def modInv(a,p):          # Finds the inverse of a mod p, if it exists
  for i in range(1,p):
    if (i*a)%p==1:
      return i
  raise ValueError(str(a)+" has no inverse mod "+str(p))

def minor(A,i,j):    # Return matrix A with the ith row and jth column deleted
  A=np.array(A)
  minor=np.zeros(shape=(len(A)-1,len(A)-1))
  p=0
  for s in range(0,len(minor)):
    if p==i:
      p=p+1
    q=0
    for t in range(0,len(minor)):
      if q==j:
        q=q+1
      minor[s][t]=A[p][q]
      q=q+1
    p=p+1
  return minor

def matrix_mod(matrix, mod_value):
    """
    Returns a new matrix where each element of 'matrix' is taken modulo 'mod_value'.
    """
    num_rows = len(matrix)

    # Initialize an empty matrix of the same dimensions

    for i in range(num_rows):
        matrix[i] = (matrix[i]) % mod_value
    return matrix


def basematrix_allon (size):
    A = [1 for _ in range(size)]
    return A
def basematrix(size):
    A = []
    print("Please enter a 1 if the light is on and a 0 if the light is off \n"
          "enter the values in order starting from the top left going right\n"
          "\t thank you!")
    for _ in range(size):
        A.append(int(input()))
    return A

def adjacency_matrix (n):
    size = n*n
    op_A = [[0 for _ in range(size)] for _ in range(size)]
    for i in range(size):
       op_A[i][i] = 1
       if size - i > n:
         op_A[i][i + n] = 1
         op_A[i + n][i] = 1
       if ((((i+1)%n)!= 0) or (i==0)) & (i < size-1):
         op_A[i+1][i] = 1
         op_A[i][i+1] = 1

    return op_A
def mod_inv(a, p):
    for x in range(1, p):
        if (a * x) % p == 1:
            return x
    raise ValueError(f"{a} has no inverse modulo {p}")

def gauss_jordan_mod_inverse(A, p):
    A = np.array(A, dtype=int)  # make sure it's an integer array
    n = A.shape[0]

    # Check A is square
    if A.shape[0] != A.shape[1]:
        raise ValueError("Matrix A must be square")

    # Augment with the identity matrix: [A | I]
    # So shape is (n, 2n)
    I = np.eye(n, dtype=int)
    #print(A%p)
    aug = np.hstack((A % p, I))  # do mod p right away
    #print(aug)
    # Perform Gauss-Jordan
    for col in range(n):
        # Find a pivot row with a non-zero pivot in 'col'
        pivot_row = col
        while pivot_row < n and aug[pivot_row, col] == 0:
            pivot_row += 1

        if pivot_row == n:
            # No pivot in this column -> matrix not invertible
            #print(aug)

            raise ValueError("Matrix is singular mod {} (no pivot in column {})".format(p, col))

        # If pivot_row != col, swap them
        if pivot_row != col:
            aug[[col, pivot_row], :] = aug[[pivot_row, col], :]


        # Now aug[col, col] is our pivot
        pivot = aug[col, col]
        inv_pivot = mod_inv(pivot, p)  # modular inverse of the pivot

        # Normalize the pivot row so pivot becomes 1
        aug[col, :] = (aug[col, :] * inv_pivot) % p

        # Eliminate this column in all other rows
        for row in range(n):
            if row != col:
                factor = aug[row, col]
                if factor != 0:
                    #print(aug, "before\n ")
                    aug[row, :] = (aug[row, :] - factor * aug[col, :]) % p
                   # print(aug,"after\n")

    # After that, left half should be the identity, and right half is A^{-1} mod p
    inverse_matrix = aug[:, n:]  # take columns n to 2n
    return inverse_matrix

def solution (n):
    mod = int(input("How many clicks on the same button does it take to return to your original state?"))
    A = basematrix(n**2)
    #A = basematrix_allon(n ** 2)

    print("here is your strategy to turn all the lights off\n"
           "click on the corresponding button where theres a 1...\n",
          matrix_mod(np.matmul(gauss_jordan_mod_inverse(adjacency_matrix(n),mod),A), mod))
          #matrix_mod(np.matmul(modMatInv(operator_matrix(n), mod), A), mod))

solution(int(input("Please enter the dimensions of the game...")))
# I've written down two version for finding the modular inverse
# Interestingly enough I can only find solutions to n dimension games where
# there are multiple solutions to the game of size n for all configurations of initial states
# some examples are 4, 5, and 9 there are multiple solutions for all so the adjacency matrix is singular

