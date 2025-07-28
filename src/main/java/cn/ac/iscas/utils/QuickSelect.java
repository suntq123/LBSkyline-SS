package cn.ac.iscas.utils;

import java.awt.font.GraphicAttribute;
import java.util.*;
import java.util.concurrent.ThreadLocalRandom;

public class QuickSelect {

    // 定义一个索引-值对的类来存储距离及其对应的索引
    static class IndexedDistance {
        double distance;
        int index;

        IndexedDistance(double distance, int index) {
            this.distance = distance;
            this.index = index;
        }
    }

    // 交换两个元素
    public static void swap(IndexedDistance[] arr, int i, int j) {
        IndexedDistance temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;
    }

    // 分区函数，返回基准元素的最终位置
    public static int partition(IndexedDistance[] arr, int left, int right, int pivotIndex) {
        double pivotValue = arr[pivotIndex].distance;
        swap(arr, pivotIndex, right);  // 将基准移到右边
        int storeIndex = left;

        // 遍历所有元素，按距离从小到大排列
        for (int i = left; i < right; i++) {
            if (arr[i].distance < pivotValue) {
                swap(arr, storeIndex, i);
                storeIndex++;
            }
        }
        swap(arr, storeIndex, right);  // 将基准放到正确的位置
        return storeIndex;
    }

    // 快速选择算法
    public static IndexedDistance[] quickSelect(IndexedDistance[] arr, int left, int right, int k) {
//        if (left > right) {
//            throw new IllegalArgumentException("Invalid range: left > right (" + left + " > " + right + ")");
//        }

        if (left == right || k == arr.length) {
            return Arrays.copyOfRange(arr, 0, k);  // 找到前 k 个最近的点
        }

        // 随机选择一个基准点
        //int pivotIndex = left + (int) (Math.random() * (right - left + 1));
        int pivotIndex = left + ThreadLocalRandom.current().nextInt(right - left + 1);
        pivotIndex = partition(arr, left, right, pivotIndex);

        // 如果基准点正好是第 k 个元素，返回前 k 个点
        if (k == pivotIndex) {
            return Arrays.copyOfRange(arr, 0, k);
        } else if (k < pivotIndex) {
            return quickSelect(arr, left, pivotIndex - 1, k);
        }
//        if (k <= pivotIndex) {
//            return quickSelect(arr, left, pivotIndex, k);
//        }
        else {
            return quickSelect(arr, pivotIndex + 1, right, k);
        }
    }

    // 计算并返回前 k 个最近的点的索引
    public static int[] findKSmallestIndexes(double[] dis, int k) {
        // 创建 IndexedDistance 数组来存储每个距离值及其对应的索引
        IndexedDistance[] arr = new IndexedDistance[dis.length];
        for (int i = 0; i < dis.length; i++) {
            arr[i] = new IndexedDistance(dis[i], i);
        }

        // 使用快速选择算法找到前 k 个最小的距离值
        IndexedDistance[] result = quickSelect(arr, 0, arr.length - 1, k);

        // 按照距离进行排序
        Arrays.sort(result, Comparator.comparingDouble(a -> a.distance));

        // 提取结果中的索引
        int[] indexes = new int[k];
        for (int i = 0; i < k; i++) {
            indexes[i] = result[i].index;
        }

        return indexes;
    }

    public static void main(String[] args) {
        // 示例数据：距离数组
        double[] dis = {10.5, 2.3, 5.8, 1.1, 9.6, 3.7, 7.4};
        int k = 7;

        // 获取前 k 个最小的距离值对应的索引
        int[] resultIndexes = findKSmallestIndexes(dis, k);

        // 输出结果
        System.out.println("前 " + k + " 个最小距离值对应的索引：");
        for (int index : resultIndexes) {
            System.out.println(index);
        }
    }
}

