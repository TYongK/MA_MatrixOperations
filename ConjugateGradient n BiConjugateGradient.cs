    public double[][] PCG(double[][] matrixA, double[][] VectorB)
    {
        int i = 0;
        int imax = matrixA.Length;
        double[][] MJacobInv = InverseJacobianMatrix(matrixA);
        double[][] xOut = MatrixCreateZero(matrixA.Length, 1);
        double[][] vectorR = MatrixSubtract(VectorB, MatrixDotProduct(matrixA, xOut));
        double[][] vectorD = MatrixDotProduct(MJacobInv, vectorR);
        double deltaNew = MatrixDotProduct(MatrixTranspose(vectorR), vectorD)[0][0];
        double deltaZero = deltaNew;

        while (i < imax && deltaNew > 0.0001*deltaZero)
        {
            double[][] Q = MatrixDotProduct(matrixA, vectorD);
            double Alpha = deltaNew / MatrixDotProduct(MatrixTranspose(vectorD), Q)[0][0];
            xOut = MatrixAddition(xOut, MatrixNumMultiply(vectorD, Alpha));
            vectorR = MatrixSubtract(vectorR,MatrixNumMultiply(Q, Alpha));
            double[][] S = MatrixDotProduct(MJacobInv, vectorR);
            double deltaOld = deltaNew;
            deltaNew = MatrixDotProduct(MatrixTranspose(vectorR), S)[0][0];
            double Beta = deltaNew / deltaOld;
            vectorD = MatrixAddition(vectorR, MatrixNumMultiply(vectorD, Beta));
            i++;
        }

        return xOut;
    }

    public double[][] PCG_Bi(double[][] matrixA, double[][] VectorB)
    {
        double[][] matrixA_Bi = MatrixCreateZero(matrixA.Length*2, matrixA.Length*2);

        for(int j = 0; j < matrixA.Length ; j++)
        {
            matrixA_Bi[j][j] = 1;
            for(int k = 0; k < matrixA.Length; k++)
            {
                matrixA_Bi[matrixA.Length + j][k] = matrixA[k][j];
                matrixA_Bi[j][matrixA.Length + k] = matrixA[j][k];
            }
        }

        double[][] VectorB_Bi = MatrixCreateZero(VectorB.Length*2, 1);
        for(int j = 0; j < VectorB.Length ; j++)
        {
            VectorB_Bi[j][0] = VectorB[j][0];
        }

        int i = 0;
        int imax = matrixA_Bi.Length;
        //double[][] MJacobInv = InverseJacobianMatrix(matrixA_Bi);
        double[][] MJacobInv = IdentityMatrix(matrixA_Bi.Length);
        double[][] xOut = MatrixCreateZero(matrixA_Bi.Length, 1);
        double[][] vectorR = MatrixSubtract(VectorB_Bi, MatrixDotProduct(matrixA_Bi, xOut));
        double[][] vectorD = MatrixDotProduct(MJacobInv, vectorR);
        double deltaNew = MatrixDotProduct(MatrixTranspose(vectorR), vectorD)[0][0];
        double deltaZero = deltaNew;


        /////deltaNew > 0.00000001*deltaZero
        while (i < imax && deltaNew > 0.00001)
        {
            double[][] Q = MatrixDotProduct(matrixA_Bi, vectorD);
            double Alpha = deltaNew / MatrixDotProduct(MatrixTranspose(vectorD), Q)[0][0];
            xOut = MatrixAddition(xOut, MatrixNumMultiply(vectorD, Alpha));
            vectorR = MatrixSubtract(vectorR,MatrixNumMultiply(Q, Alpha));
            double[][] S = MatrixDotProduct(MJacobInv, vectorR);
            double deltaOld = deltaNew;
            deltaNew = MatrixDotProduct(MatrixTranspose(vectorR), S)[0][0];
            double Beta = deltaNew / deltaOld;
            vectorD = MatrixAddition(vectorR, MatrixNumMultiply(vectorD, Beta));
            i++;
        }



        double[][] xOut_Bi = MatrixCreateZero(xOut.Length/2, 1);

        for(int j = 0; j < xOut_Bi.Length ; j++)
        {
            xOut_Bi[j][0] = xOut[xOut_Bi.Length+j][0];
        }


        return xOut_Bi;
    }

    static double[][] MatrixCreate(int rows, int cols)
    {
        double[][] result = new double[rows][];

        for (int i = 0; i < rows; ++i)
        {
        result[i] = new double[cols];
        }
        return result;
    }

    static double[][] IdentityMatrix(int rowNCol)
    {
        double[][] result = new double[rowNCol][];

        for (int i = 0; i < rowNCol; ++i)
        {
            result[i] = new double[rowNCol];
        }

        for (int i = 0; i < rowNCol; i++)
        {
            for (int j = 0; j < rowNCol; j++)
            {
                if(i == j)
                {
                    result[i][j] = 1;
                }
                else
                {
                    result[i][j] = 0;
                }
            }
        }

        return result;
    }

    static double[][] InverseJacobianMatrix(double[][] matrix)
    {
        int rowNCol = matrix.Length;
        double[][] result = new double[rowNCol][];

        for (int i = 0; i < rowNCol; ++i)
        {
            result[i] = new double[rowNCol];
        }

        for (int i = 0; i < rowNCol; i++)
        {
            for (int j = 0; j < rowNCol; j++)
            {
                if(i == j)
                {
                    result[i][j] = 1 / matrix[i][j];
                }
                else
                {
                    result[i][j] = 0;
                }
            }
        }

        return result;
    }

    static double[][] MatrixCreateZero(int rows, int cols)
    {
        double[][] result = new double[rows][];

        for (int i = 0; i < rows; ++i)
        {
            result[i] = new double[cols];
        }

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                result[i][j] = 0;
            }
        }

        return result;
    }

    static double[][] MatrixDuplicate(double[][] matrix)
    {
        // allocates/creates a duplicate of a matrix.
        double[][] result = MatrixCreate(matrix.Length, matrix[0].Length);
        for (int i = 0; i < matrix.Length; ++i) // copy the values
        {
        for (int j = 0; j < matrix[i].Length; ++j)
        {
            result[i][j] = matrix[i][j];
        }
        }
        return result;
    }

    static double[][] MatrixTranspose(double[][] matrix)
    {
        // allocates/creates a duplicate of a matrix.
        double[][] result = MatrixCreate(matrix[0].Length, matrix.Length);
        for (int i = 0; i < matrix.Length; ++i) // copy the values
        {
            for (int j = 0; j < matrix[i].Length; ++j)
            {
                result[j][i] = matrix[i][j];
            }
        }
        return result;
    }

    static dynamic MatrixDotProduct(double[][] matrixA, double[][] matrixB)
    {
        int aRows = matrixA.Length; int aCols = matrixA[0].Length;
        int bRows = matrixB.Length; int bCols = matrixB[0].Length;
        if (aCols != bRows)
        {
        return "Non-conformable matrices";
        }

        double[][] result = MatrixCreate(aRows, bCols);

        Parallel.For(0, aRows, i =>
        {
        //for (int i = 0; i < aRows; ++i) // each row of A
            for (int j = 0; j < bCols; ++j) // each col of B
                for (int k = 0; k < aCols; ++k) // could use k less-than bRows
                result[i][j] += matrixA[i][k] * matrixB[k][j];
        });
        return result;
    }

    static dynamic MatrixSubtract(double[][] matrixA, double[][] matrixB)
    {
        int aRows = matrixA.Length; int aCols = matrixA[0].Length;
        int bRows = matrixB.Length; int bCols = matrixB[0].Length;
        if (aRows != bRows || aCols != bCols)
        {
            return "Non-conformable matrices";
        }

        double[][] result = MatrixCreate(aRows, bCols);

        //for (int i = 0; i < aRows; ++i)
        Parallel.For(0, aRows, i =>
        {
            for (int j = 0; j < bCols; ++j)
                result[i][j] = matrixA[i][j] - matrixB[i][j];
        });
        return result;
    }

    static dynamic MatrixAddition(double[][] matrixA, double[][] matrixB)
    {
        int aRows = matrixA.Length; int aCols = matrixA[0].Length;
        int bRows = matrixB.Length; int bCols = matrixB[0].Length;
        if (aRows != bRows || aCols != bCols)
        {
            return "Non-conformable matrices";
        }

        double[][] result = MatrixCreate(aRows, bCols);

        //for (int i = 0; i < aRows; ++i)
        Parallel.For(0, aRows, i =>
        {
        for (int j = 0; j < bCols; ++j)
            result[i][j] = matrixA[i][j] + matrixB[i][j];
        });

        return result;
    }

    static dynamic MatrixNumMultiply(double[][] matrixA, double numIn)
    {
        int aRows = matrixA.Length; int aCols = matrixA[0].Length;

        double[][] result = MatrixCreate(aRows, aCols);

        //for (int i = 0; i < aRows; ++i)
        Parallel.For(0, aRows, i =>
        {
        for (int j = 0; j < aCols; ++j)
            result[i][j] = matrixA[i][j] * numIn;
        });

        return result;
    }
