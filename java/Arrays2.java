/**
 * Arrays2
 *
 * @author Ludo Waltman
 * @author Nees Jan van Eck
 * @version 1.3.1, 11/17/14
 */

import java.util.Arrays;
import java.util.Random;

public class Arrays2
{
    public static double calcSum(double[] value)
    {
        double sum;
        int i;

        sum = 0;
        for (i = 0; i < value.length; i++)
            sum += value[i];
        return sum;
    }

    public static double calcSum(double[] value, int beginIndex, int endIndex)
    {
        double sum;
        int i;

        sum = 0;
        for (i = beginIndex; i < endIndex; i++)
            sum += value[i];
        return sum;
    }

    public static double calcAverage(double[] value)
    {
        double average;
        int i;

        average = 0;
        for (i = 0; i < value.length; i++)
            average += value[i];
        average /= value.length;
        return average;
    }

    public static double calcMedian(double[] value)
    {
        double median;
        double[] sortedValue;

        sortedValue = (double[])value.clone();
        Arrays.sort(sortedValue);
        if (sortedValue.length % 2 == 1)
            median = sortedValue[(sortedValue.length - 1) / 2];
        else
            median = (sortedValue[sortedValue.length / 2 - 1] + sortedValue[sortedValue.length / 2]) / 2;
        return median;
    }

    public static double calcMinimum(double[] value)
    {
        double minimum;
        int i;

        minimum = value[0];
        for (i = 1; i < value.length; i++)
            minimum = Math.min(minimum, value[i]);
        return minimum;
    }

    public static double calcMaximum(double[] value)
    {
        double maximum;
        int i;

        maximum = value[0];
        for (i = 1; i < value.length; i++)
            maximum = Math.max(maximum, value[i]);
        return maximum;
    }

    public static int calcMaximum(int[] value)
    {
        int i, maximum;

        maximum = value[0];
        for (i = 1; i < value.length; i++)
            maximum = Math.max(maximum, value[i]);
        return maximum;
    }

    public static int[] generateRandomPermutation(int nElements)
    {
        return generateRandomPermutation(nElements, new Random());
    }

    public static int[] generateRandomPermutation(int nElements, Random random)
    {
        int i, j, k;
        int[] permutation;

        permutation = new int[nElements];
        for (i = 0; i < nElements; i++)
            permutation[i] = i;
        for (i = 0; i < nElements; i++)
        {
            j = random.nextInt(nElements);
            k = permutation[i];
            permutation[i] = permutation[j];
            permutation[j] = k;
        }
        return permutation;
    }
}
