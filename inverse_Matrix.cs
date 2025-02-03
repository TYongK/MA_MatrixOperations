using System;
using System.Collections;
using System.Collections.Generic;

using Rhino;
using Rhino.Geometry;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;



/// <summary>
/// This class will be instantiated on demand by the Script component.
/// </summary>
public class Script_Instance : GH_ScriptInstance
{
#region Utility functions
  /// <summary>Print a String to the [Out] Parameter of the Script component.</summary>
  /// <param name="text">String to print.</param>
  private void Print(string text) { /* Implementation hidden. */ }
  /// <summary>Print a formatted String to the [Out] Parameter of the Script component.</summary>
  /// <param name="format">String format.</param>
  /// <param name="args">Formatting parameters.</param>
  private void Print(string format, params object[] args) { /* Implementation hidden. */ }
  /// <summary>Print useful information about an object instance to the [Out] Parameter of the Script component. </summary>
  /// <param name="obj">Object instance to parse.</param>
  private void Reflect(object obj) { /* Implementation hidden. */ }
  /// <summary>Print the signatures of all the overloads of a specific method to the [Out] Parameter of the Script component. </summary>
  /// <param name="obj">Object instance to parse.</param>
  private void Reflect(object obj, string method_name) { /* Implementation hidden. */ }
#endregion

#region Members
  /// <summary>Gets the current Rhino document.</summary>
  private readonly RhinoDoc RhinoDocument;
  /// <summary>Gets the Grasshopper document that owns this script.</summary>
  private readonly GH_Document GrasshopperDocument;
  /// <summary>Gets the Grasshopper script component that owns this script.</summary>
  private readonly IGH_Component Component;
  /// <summary>
  /// Gets the current iteration count. The first call to RunScript() is associated with Iteration==0.
  /// Any subsequent call within the same solution will increment the Iteration count.
  /// </summary>
  private readonly int Iteration;
#endregion

  /// <summary>
  /// This procedure contains the user code. Input parameters are provided as regular arguments,
  /// Output parameters as ref arguments. You don't have to assign output parameters,
  /// they will have a default value.
  /// </summary>
  private void RunScript(DataTree<double> Tree_axa_Matrix, ref object A)
  {

    if(MatrixInverse_Grasshopper(Tree_axa_Matrix).GetType() == typeof(string))
    {
      A = "Non-conformable matrices";
    }
    else
    {
      double[][] m = MatrixInverse_Grasshopper(Tree_axa_Matrix);

      //////Define Tree//////
      DataTree < double > doubleTree = new DataTree<double>();

      // Raws = Branches
      for(int i = 0; i < m.Length;i++)
      {
        GH_Path pth = new GH_Path(i);

        // Cols = Index
        for (int j = 0; j < m[i].Length;j++)
        {
          doubleTree.Add(m[i][j], pth);
        }
      }
      //////Define Tree End//////

      A = doubleTree;
    }


  }

  // <Custom additional code> 
  //////////////////////////////////
  static double[][] MatrixCreate(int rows, int cols)
  {
    double[][] result = new double[rows][];

    for (int i = 0; i < rows; ++i)
    {
      result[i] = new double[cols];
    }
    return result;
  }
  //////////////////////////////////

  //////////////////////////////////
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
  //////////////////////////////////

  //////////////////////////////////
  static double[][] MatrixDecompose_Grasshopper(DataTree<double> matrix, out int[] perm, out int toggle, out string status)
  {
    // Doolittle LUP decomposition with partial pivoting.
    // rerturns: result is L (with 1s on diagonal) and U;
    // perm holds row permutations; toggle is +1 or -1 (even or odd)
    IList<List<double>> matrix_Branches = matrix.Branches;
    int rows = matrix_Branches.Count; int cols = matrix_Branches[0].Count; // assume square

    status = "Works fine";

    if (rows != cols)
    {
      status = "Attempt to decompose a non-square m";
    }


    //Tree into [][] Start//
    double[][] matrix_TtoM = MatrixCreate(rows, cols);
    for(int i = 0; i < rows; i++)
    {
      for(int j = 0; j < cols; j++)
      {
        matrix_TtoM[i][j] = matrix_Branches[i][j];
      }
    }
    //Tree into [][] End//

    int n = rows; // convenience

    double[][] result = MatrixDuplicate(matrix_TtoM);

    perm = new int[n]; // set up row permutation result

    for (int i = 0; i < n; ++i)
    {
      perm[i] = i;
    }

    toggle = 1; // toggle tracks row swaps.
    // +1 -greater-than even, -1 -greater-than odd. used by MatrixDeterminant

    for (int j = 0; j < n - 1; ++j) // each column
    {
      double colMax = Math.Abs(result[j][j]); // find largest val in col
      int pRow = j;
      //for (int i = j + 1; i less-than n; ++i)
      //{
      //  if (result[i][j] greater-than colMax)
      //  {
      //    colMax = result[i][j];
      //    pRow = i;
      //  }
      //}

      // reader Matt V needed this:
      for (int i = j + 1; i < n; ++i)
      {
        if (Math.Abs(result[i][j]) > colMax)
        {
          colMax = Math.Abs(result[i][j]);
          pRow = i;
        }
      }
      // Not sure if this approach is needed always, or not.

      if (pRow != j) // if largest value not on pivot, swap rows
      {
        double[] rowPtr = result[pRow];
        result[pRow] = result[j];
        result[j] = rowPtr;

        int tmp = perm[pRow]; // and swap perm info
        perm[pRow] = perm[j];
        perm[j] = tmp;

        toggle = -toggle; // adjust the row-swap toggle
      }

      // --------------------------------------------------
      // This part added later (not in original)
      // and replaces the 'return null' below.
      // if there is a 0 on the diagonal, find a good row
      // from i = j+1 down that doesn't have
      // a 0 in column j, and swap that good row with row j
      // --------------------------------------------------

      if (result[j][j] == 0.0)
      {
        // find a good row to swap
        int goodRow = -1;
        for (int row = j + 1; row < n; ++row)
        {
          if (result[row][j] != 0.0)
            goodRow = row;
        }

        if (goodRow == -1)
        {
          status = "Cannot use Doolittle's method";
        }

        // swap rows so 0.0 no longer on diagonal
        double[] rowPtr = result[goodRow];
        result[goodRow] = result[j];
        result[j] = rowPtr;

        int tmp = perm[goodRow]; // and swap perm info
        perm[goodRow] = perm[j];
        perm[j] = tmp;

        toggle = -toggle; // adjust the row-swap toggle
      }
      // --------------------------------------------------
      // if diagonal after swap is zero . .
      //if (Math.Abs(result[j][j]) less-than 1.0E-20)
      //  return null; // consider a throw

      for (int i = j + 1; i < n; ++i)
      {
        result[i][j] /= result[j][j];
        for (int k = j + 1; k < n; ++k)
        {
          result[i][k] -= result[i][j] * result[j][k];
        }
      }


    } // main j column loop

    return result;
  } // MatrixDecompose
  //////////////////////////////////

  //////////////////////////////////
  static double[] HelperSolve(double[][] luMatrix, double[] b)
  {
    // before calling this helper, permute b using the perm array
    // from MatrixDecompose that generated luMatrix
    int n = luMatrix.Length;
    double[] x = new double[n];
    b.CopyTo(x, 0);

    for (int i = 1; i < n; ++i)
    {
      double sum = x[i];
      for (int j = 0; j < i; ++j)
      {
        sum -= luMatrix[i][j] * x[j];
      }
      x[i] = sum;
    }

    x[n - 1] /= luMatrix[n - 1][n - 1];
    for (int i = n - 2; i >= 0; --i)
    {
      double sum = x[i];
      for (int j = i + 1; j < n; ++j)
      {
        sum -= luMatrix[i][j] * x[j];
      }
      x[i] = sum / luMatrix[i][i];
    }

    return x;
  }
  //////////////////////////////////

  //////////////////////////////////
  //returns "string" or "double[][]"
  static dynamic MatrixInverse_Grasshopper(DataTree<double> matrix)
  {
    IList<List<double>> matrix_Branches = matrix.Branches;
    int rows = matrix_Branches.Count; int cols = matrix_Branches[0].Count; // assume square

    //Tree into [][] Start//
    double[][] matrix_TtoM = MatrixCreate(rows, cols);
    for(int i = 0; i < rows; i++)
    {
      for(int j = 0; j < cols; j++)
      {
        matrix_TtoM[i][j] = matrix_Branches[i][j];
      }
    }
    //Tree into [][] End//

    int n = matrix.BranchCount;
    double[][] result = matrix_TtoM;

    int[] perm;
    int toggle;
    string state;
    double[][] lum = MatrixDecompose_Grasshopper(matrix, out perm, out toggle, out state);

    if (lum == null)
    {
      return "Unable to compute inverse";
    }

    double[] b = new double[n];
    for (int i = 0; i < n; ++i)
    {
      for (int j = 0; j < n; ++j)
      {
        if (i == perm[j])
          b[j] = 1.0;
        else
          b[j] = 0.0;
      }

      double[] x = HelperSolve(lum, b); //

      for (int j = 0; j < n; ++j)
        result[j][i] = x[j];
    }
    return result;
  }

  static dynamic MatrixInverse(double[][] matrix)
  {
      int rows = matrix.Length;
      int cols = matrix[0].Length; // assume square

      //Tree into [][] Start//
      double[][] matrix_TtoM = MatrixCreate(rows, cols);
      for(int i = 0; i < rows; i++)
      {
          for(int j = 0; j < cols; j++)
          {
              matrix_TtoM[i][j] = matrix[i][j];
          }
      }
      //Tree into [][] End//

      int n = matrix.Length;
      double[][] result = matrix_TtoM;

      int[] perm;
      int toggle;
      string state;
      double[][] lum = MatrixDecompose(matrix, out perm, out toggle, out state);

      if (lum == null)
      {
      return "Unable to compute inverse";
      }

      double[] b = new double[n];
      for (int i = 0; i < n; ++i)
      {
          for (int j = 0; j < n; ++j)
          {
              if (i == perm[j])
              b[j] = 1.0;
              else
              b[j] = 0.0;
          }

          double[] x = HelperSolve(lum, b); //

          for (int j = 0; j < n; ++j)
              result[j][i] = x[j];
      }
      return result;
  }
  //////////////////////////////////
  // Resource#1: https://jamesmccaffrey.wordpress.com/2015/03/06/inverting-a-matrix-using-c/
  // Resource#2: https://en.wikipedia.org/wiki/LU_decomposition
  // Resource#3: https://www.youtube.com/watch?v=dza5JTvMpzk
  // Resource#4: http://home.cc.umanitoba.ca/~farhadi/Math2120/Inverse%20Using%20LU%20decomposition.pdf
  // Editor: Taeyong Kim.
  // </Custom additional code> 
}
