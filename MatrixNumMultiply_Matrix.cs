    static dynamic MatrixNumMultiply(double[][] matrixA, double numIn)
    {
        int aRows = matrixA.Length; int aCols = matrixA[0].Length;

        double[][] result = MatrixCreate(aRows, aCols);

        for (int i = 0; i < aRows; ++i)
        for (int j = 0; j < aCols; ++j)
            result[i][j] = matrixA[i][j] * numIn;

        return result;
    }
