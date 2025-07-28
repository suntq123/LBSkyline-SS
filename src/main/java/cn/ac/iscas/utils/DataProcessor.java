package cn.ac.iscas.utils;

import cn.ac.iscas.kdtree.KDTreeNode;

import java.io.*;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.*;
import java.util.AbstractMap.SimpleEntry;

import cn.ac.iscas.skylinequery.SkylineQuery;
import com.alibaba.fastjson.JSONObject;

import static cn.ac.iscas.secretsharing.AdditiveSecretSharing.multiply;
import static cn.ac.iscas.secretsharing.AdditiveSecretSharing.randomSplit;

public class DataProcessor {

    public static BigInteger[][] generateDataset(int number, int dimension, int bitLength, Random random) {
        BigInteger[][] dataset = new BigInteger[number][dimension + 1];

        // Random random = new Random();
        for (int i = 0; i < number; i++) {
            for (int j = 0; j < dimension; j++) {
                dataset[i][j] = new BigInteger(bitLength, random);
            }
            dataset[i][dimension] = BigInteger.valueOf(i + 1); // ptr
        }

        return dataset;
    }

    public static BigInteger[][] generateQueryPoints_HOTE(int number, int dimension, int bitLength, Random random) {
        BigInteger[][] dataset = new BigInteger[number][dimension + 1];

        // Random random = new Random();
        for (int i = 0; i < number; i++) {
            for (int j = 0; j < dimension; j++) {
                dataset[i][j] = new BigInteger(bitLength, random);
            }
            dataset[i][dimension] = BigInteger.valueOf(i + 1); // ptr
        }

        return dataset;
    }

    public static BigInteger[][] generateDataset_CORR(int number, int dimension, int bitLength, Random random) {
        BigInteger[][] dataset = new BigInteger[number][dimension + 1];
        double coe = 0.8;

        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < number; j++) {
                if(i<3) {
                    dataset[j][i] = new BigInteger(bitLength, random);
                }else {
                    Random randomi = new Random(1000L + j + (long) (i-1) *number);
                    dataset[j][i] = dataset[j][i - 1].add(BigInteger.valueOf((long)(coe * randomi.nextDouble()*10))).mod(BigInteger.TWO.pow(bitLength));
                }
            }
        }
        for (int j = 0; j < number; j++) {
            dataset[j][dimension] = BigInteger.valueOf(j + 1);
        }

        return dataset;
    }

    public static BigInteger[][] generateDataset_ANTI(int number, int dimension, int bitLength, Random random) {
        BigInteger[][] dataset = new BigInteger[number][dimension + 1];
        double coe = 0.8;

        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < number; j++) {
                if(i<3) {
                    dataset[j][i] = new BigInteger(bitLength, random);
                }else {
                    Random randomi = new Random(1000L + j + (long) (i-1) *number);
                    dataset[j][i] = BigInteger.TWO.pow(bitLength).subtract(BigInteger.valueOf((long)(dataset[j][i - 1].intValue()*coe))).add(BigInteger.valueOf((long)(coe * randomi.nextDouble()*10))).mod(BigInteger.TWO.pow(bitLength));
                }
            }
        }
        for (int j = 0; j < number; j++) {
            dataset[j][dimension] = BigInteger.valueOf(j + 1);
        }

        return dataset;
    }

    public static BigInteger[][] generateDataset_INDE(int number, int dimension, int bitLength, Random random) {
        BigInteger[][] dataset = new BigInteger[number][dimension + 1];

        for (int i = 0; i < number; i++) {
            for (int j = 0; j < dimension; j++) {
                dataset[i][j] = new BigInteger(bitLength, random);
            }
            dataset[i][dimension] = BigInteger.valueOf(i + 1); // ptr
        }

        return dataset;
    }

    public static BigInteger[][] getDataset_HOTE(int number, int dimension, int bitLength, Random random) {
        BigInteger[][] dataset = new BigInteger[number][dimension + 1];

        String[] mmtStr = Util.readTxtFile("src\\main\\data\\makemytrip_com-travel_sample.txt");
        for (int i = 0; i < number; i++) {
            String[] hotelStr = mmtStr[i].trim().split(" ");
            for (int j = 0; j < dimension; j++) {
                dataset[i][j] = new BigInteger(hotelStr[j]);
            }
            dataset[i][dimension] = new BigInteger(hotelStr[7]);
        }

        return dataset;
    }

    public static BigInteger[][] generateQueryDataset(int number, int dimension, int bitLength, double area, Random random) {
        BigInteger[][] dataset = new BigInteger[number][dimension + 1];

        for (int j = 0; j < dimension; j++) {
            dataset[0][j] = new BigInteger(bitLength, random);
        }
        dataset[0][dimension] = BigInteger.valueOf(0); // ptr
        int bitLengthV = (int) (Util.log2(Math.sqrt(Math.pow(2, bitLength) * Math.pow(2, bitLength) * area)));
        for (int i = 1; i < number; i++) {
            for (int j = 0; j < dimension; j++) {
                dataset[i][j] = (dataset[0][j].add(new BigInteger(bitLengthV, random))).mod(BigInteger.TWO.pow(bitLength));
            }
            dataset[i][dimension] = BigInteger.valueOf(i + 1); // ptr
        }

        return dataset;
    }

    public static BigInteger calDist(BigInteger[] x, BigInteger[] y) {

        BigInteger sum = BigInteger.ZERO;
        for (int i = 0; i < 2; i++) { //x.length
            sum = sum.add((x[i].subtract(y[i])).multiply(x[i].subtract(y[i])));
        }

        return sum.sqrt();
    }

    public static BigInteger sumArray(BigInteger[] point) {
        BigInteger sum = BigInteger.ZERO;
        for (int i = 0; i < point.length; i++) {
            sum = sum.add(point[i]);
        }
        return sum;
    }

    public static int min(List<BigInteger> list) {
        BigInteger min = list.get(0);
        int index = 0;
        for (int i = 1; i < list.size(); i++) {
            if (list.get(i).compareTo(min) < 0) {
                min = list.get(i);
                index = i;
            }
        }
        return index;
    }

    public static boolean dominanceCheck(BigInteger[] x, BigInteger[] y) {
        int countEquals = 0;
        for (int i = 0; i < x.length; i++) {

            if (x[i].compareTo(y[i]) > 0) {
                return false;
            }else if (x[i].compareTo(y[i]) == 0) {
                    countEquals++;
                }
        }
        if (countEquals == x.length) {
            return false;
        }
        return true;
    }

    public static boolean dominanceCheck2(BigInteger[] x, BigInteger[] y) {
        for (int i = 0; i < x.length; i++) {

            if (x[i].compareTo(y[i]) > 0) {
                return false;
            }
        }
        return true;
    }

    public static Set<SkylineQuery.Point> getSkyline(BigInteger[][] dataset, BigInteger[][] query,
                                                     int dataNumber, int dimension, int q){
        List<BigInteger[]> queryList = new ArrayList<>();
        List<BigInteger[]> pointList = new ArrayList<>();
        for (int i = 0; i < dataNumber; i++) {
//            BigInteger[] point = new BigInteger[dimension + 1];
//            for (int j = 0; j < dimension + 1; j++) {
//                point[j] = dataset[i][j];
//            }
//            pointList.add(point);
            pointList.add(dataset[i]);
        }
        for (int i = 0; i < q; i++) {
            BigInteger[] point = new BigInteger[2];
            for (int j = 0; j < 2; j++) {
                point[j] = query[i][j];
            }
            queryList.add(point);
        }

        List<BigInteger[]> skylineList = new ArrayList<>();
        List<BigInteger[]> tempResult1 = new ArrayList<>();
        List<BigInteger> tempResult2 = new ArrayList<>();

        if (queryList != null && pointList != null)
        {
            for (int i = 0; i < pointList.size(); i++)
            {
                BigInteger[] point = pointList.get(i);
                BigInteger[] temp = new BigInteger[point.length-1-2+queryList.size()];

                for(int k=0; k<queryList.size(); k++){
                    temp[k] = calDist(queryList.get(k), point);
                }

                for (int j = queryList.size(); j < temp.length; j++)
                {
                    temp[j] = point[j-(queryList.size()-2)];
                }
                tempResult1.add(temp);
            }
            for (BigInteger[] point : tempResult1) {
                BigInteger temp = sumArray(point);
                tempResult2.add(temp);
            }
            while (!tempResult2.isEmpty()){
                int res= min(tempResult2);
                BigInteger[] skylinePoint = pointList.get(res);

//                for (int i = 0; i < skylinePoint.length; i++){
//                    System.out.print(skylinePoint[i]+" ");
//                }
//                System.out.println();

                skylineList.add(skylinePoint);
                BigInteger[] skylinePointPrime = tempResult1.get(res);
                pointList.remove(res);
                tempResult1.remove(res);
                tempResult2.remove(res);

                List<BigInteger[]> delpoints = new ArrayList<>();
                for(BigInteger[] point : tempResult1){
                    if(dominanceCheck(skylinePointPrime,point)){
                        delpoints.add(point);
                    }
                }
                for(int i = 0; i < delpoints.size(); i++){
                    BigInteger[] point = delpoints.get(i);
                    int index = tempResult1.indexOf(point);
                    pointList.remove(index);
                    tempResult1.remove(index);
                    tempResult2.remove(index);
                }
            }

        }
        Set<SkylineQuery.Point> result = new HashSet<>();
        for(BigInteger[] skylinePoint : skylineList){
            SkylineQuery.Point point = new SkylineQuery.Point(dimension);
            for(int i = 0; i < dimension; i++){
                point.data[i] = skylinePoint[i];
            }
            point.id = skylinePoint[dimension];
            result.add(point);
        }
        return result;
    }

    public static Set<List> getNonspatialDominators(BigInteger[][] dataset, int dataNumber, int dimension){

        List<SkylineQuery.Point> pointList = new ArrayList<>();
        for (int i = 0; i < dataNumber; i++) {
            SkylineQuery.Point point = new SkylineQuery.Point(dimension);
            point.id = dataset[i][dimension];
            for (int j = 0; j < dimension; j++) {
                point.data[j] = dataset[i][j];
            }
            pointList.add(point);
        }

        Set<List> result = new HashSet<>();

        for (int i = 0; i < dataNumber; i++)
        {
            SkylineQuery.Point pointi = pointList.get(i);
            List<SkylineQuery.Point> dominatorsList = new ArrayList<>();
            dominatorsList.add(pointi);
            for (int j = 0; j < dataNumber; j++)
            {
                SkylineQuery.Point pointj = pointList.get(j);
                if(dominanceCheck(Arrays.copyOfRange(pointj.data,2,dimension),Arrays.copyOfRange(pointi.data,2,dimension))){
                    dominatorsList.add(pointj);
                }
            }
            result.add(dominatorsList);
        }

        return result;
    }

    public static List<SkylineQuery.Point> parseDataset(BigInteger[][] dataset, int dimension){
        List<SkylineQuery.Point> pointList = new ArrayList<>();
        for (int i = 0; i < dataset.length; i++) {
            SkylineQuery.Point point = new SkylineQuery.Point(dimension);
            point.id = dataset[i][dimension];
            for (int j = 0; j < dimension; j++) {
                point.data[j] = dataset[i][j];
            }
            pointList.add(point);
        }
        return pointList;
    }

    public static BigInteger[][] pointToDataset(Set<SkylineQuery.Point> pointset, int dimension){
        BigInteger[][] dataset = new BigInteger[pointset.size()][dimension+1];
        int index = 0;
        for (SkylineQuery.Point point : pointset){
            dataset[index][dimension] = point.id;
            for (int j = 0; j < dimension; j++){
                dataset[index][j] = point.data[j];
            }
            index++;
        }
        return dataset;
    }

    public static List<SkylineQuery.Point> getNonspatialDominators2(List<SkylineQuery.Point> pointList, SkylineQuery.Point point, int id, int dimension){

        pointList.remove(id);
        List<SkylineQuery.Point> dominatorsList = new ArrayList<>();
//        dominatorsList.add(point);

        for (int i = 0; i < pointList.size(); i++)
        {
            SkylineQuery.Point pointi = pointList.get(i);
            if(dominanceCheck2(Arrays.copyOfRange(pointi.data,2,dimension),Arrays.copyOfRange(point.data,2,dimension))){
                dominatorsList.add(pointi);
            }
        }

        return dominatorsList;
    }

    public static List<SkylineQuery.Point> getNonspatialDominators3(List<SkylineQuery.Point> pointList, SkylineQuery.Point point, int dimension){

        List<SkylineQuery.Point> dominatorsList = new ArrayList<>();
        List<SkylineQuery.Point> dominatorsListR = new ArrayList<>();
        List<BigInteger> distList = new ArrayList<>();
        //dominatorsList.add(point);

        for (int i = 0; i < pointList.size(); i++)
        {
            SkylineQuery.Point pointi = pointList.get(i);
            if(dominanceCheck2(Arrays.copyOfRange(pointi.data,2,dimension),Arrays.copyOfRange(point.data,2,dimension))){
                dominatorsList.add(pointi);
                BigInteger dist  = calDist(Arrays.copyOfRange(pointi.data,0,2),Arrays.copyOfRange(point.data,0,2));
                distList.add(dist);
            }
        }

        List<Integer> indices = new ArrayList<>();
        for (int i = 0; i < distList.size(); i++) {
            indices.add(i);
        }

        indices.sort(Comparator.comparing(i -> distList.get(i)));
        for (int i = 0; i < indices.size(); i++) {
            int index = indices.get(i);
            dominatorsListR.add(dominatorsList.get(index));
        }

        return dominatorsListR;
    }

    public static boolean isIntersect(List<SkylineQuery.Vertex> polygon1, SkylineQuery.Point low, SkylineQuery.Point high) {
        boolean isIntersect = false;
        for (int i = 0; i < polygon1.size(); i++) {
            SkylineQuery.Vertex vertex1 = polygon1.get(i);
            SkylineQuery.Vertex vertex2 = polygon1.get((i + 1) % polygon1.size());

            if (vertex1.data[0].compareTo(vertex2.data[0]) == 0 && vertex1.data[1].compareTo(vertex2.data[1]) == 0) {
                continue;
            }
            if (((high.data[0].subtract(low.data[0]).multiply(vertex1.data[1].subtract(high.data[1]))).multiply(high.data[0].subtract(low.data[0]).multiply(vertex2.data[1].subtract(high.data[1]))).compareTo(BigInteger.ZERO) < 0
            && ((vertex2.data[0].subtract(vertex1.data[0]).multiply(high.data[1].subtract(vertex1.data[1]))).subtract(vertex2.data[1].subtract(vertex1.data[1]).multiply(low.data[0].subtract(vertex1.data[0])))).multiply(
               ((vertex2.data[0].subtract(vertex1.data[0]).multiply(high.data[1].subtract(vertex1.data[1]))).subtract(vertex2.data[1].subtract(vertex1.data[1]).multiply(high.data[0].subtract(vertex1.data[0]))))).compareTo(BigInteger.ZERO) <= 0)
            || ((low.data[1].subtract(high.data[1]).multiply(vertex1.data[0].subtract(high.data[0]))).multiply(low.data[1].subtract(high.data[1]).multiply(vertex2.data[0].subtract(high.data[0]))).compareTo(BigInteger.ZERO) < 0
            && ((vertex2.data[0].subtract(vertex1.data[0]).multiply(high.data[1].subtract(vertex1.data[1]))).subtract(vertex2.data[1].subtract(vertex1.data[1]).multiply(high.data[0].subtract(vertex1.data[0])))).multiply(
               ((vertex2.data[0].subtract(vertex1.data[0]).multiply(low.data[1].subtract(vertex1.data[1]))).subtract(vertex2.data[1].subtract(vertex1.data[1]).multiply(high.data[0].subtract(vertex1.data[0]))))).compareTo(BigInteger.ZERO) <= 0)
            || ((low.data[0].subtract(high.data[0]).multiply(vertex1.data[1].subtract(low.data[1]))).multiply(low.data[0].subtract(high.data[0]).multiply(vertex2.data[1].subtract(low.data[1]))).compareTo(BigInteger.ZERO) < 0
            && ((vertex2.data[0].subtract(vertex1.data[0]).multiply(low.data[1].subtract(vertex1.data[1]))).subtract(vertex2.data[1].subtract(vertex1.data[1]).multiply(high.data[0].subtract(vertex1.data[0])))).multiply(
               ((vertex2.data[0].subtract(vertex1.data[0]).multiply(low.data[1].subtract(vertex1.data[1]))).subtract(vertex2.data[1].subtract(vertex1.data[1]).multiply(low.data[0].subtract(vertex1.data[0]))))).compareTo(BigInteger.ZERO) <= 0)
            || ((high.data[1].subtract(low.data[1]).multiply(vertex1.data[0].subtract(low.data[0]))).multiply(high.data[1].subtract(low.data[1]).multiply(vertex2.data[0].subtract(low.data[0]))).compareTo(BigInteger.ZERO) < 0
            && ((vertex2.data[0].subtract(vertex1.data[0]).multiply(low.data[1].subtract(vertex1.data[1]))).subtract(vertex2.data[1].subtract(vertex1.data[1]).multiply(low.data[0].subtract(vertex1.data[0])))).multiply(
               ((vertex2.data[0].subtract(vertex1.data[0]).multiply(high.data[1].subtract(vertex1.data[1]))).subtract(vertex2.data[1].subtract(vertex1.data[1]).multiply(low.data[0].subtract(vertex1.data[0]))))).compareTo(BigInteger.ZERO) <= 0)){
                isIntersect = true;
            }
        }
        return isIntersect;
    }

    public static Set<BigInteger> getKNearest(BigInteger[][] dataset, BigInteger[] point,
            int dimension, int pow, int k) {
        Set<BigInteger> result = new HashSet<>();

        // <pId, distance>
        List<SimpleEntry<BigInteger, BigInteger>> knnList = new ArrayList<>();
        for (int i = 0; i < dataset.length; i++) {
            BigInteger id = dataset[i][dimension];

            BigInteger d = BigInteger.ZERO;
            for (int j = 0; j < dimension; j++) {
                BigInteger absDis = dataset[i][j].subtract(point[j]).abs();
                if (pow <= 0) {
                    if (d.compareTo(absDis) < 0)
                        d = absDis;
                } else {
                    d = d.add(absDis.pow(pow));
                }
            }

            knnList.add(new SimpleEntry<BigInteger, BigInteger>(id, d));
        }

        knnList.sort((e1, e2) -> {
            return e1.getValue().compareTo(e2.getValue());
        });

        for (int i = 0; i < k; i++) {
            result.add(knnList.get(i).getKey());
        }

        return result;
    }

    public static void main(String[] args) throws IOException {

        SkylineQuery.Point low = new SkylineQuery.Point(BigInteger.valueOf(0), new BigInteger[] { BigInteger.valueOf(0), BigInteger.valueOf(0) });
        SkylineQuery.Point high = new SkylineQuery.Point(BigInteger.valueOf(0), new BigInteger[] { BigInteger.valueOf(4), BigInteger.valueOf(4) });
        List<SkylineQuery.Vertex> vertices = new ArrayList<>();
        vertices.add(new SkylineQuery.Vertex(new BigInteger[] { BigInteger.valueOf(1), BigInteger.valueOf(3) }));
        vertices.add(new SkylineQuery.Vertex(new BigInteger[] { BigInteger.valueOf(2), BigInteger.valueOf(5) }));
        vertices.add(new SkylineQuery.Vertex(new BigInteger[] { BigInteger.valueOf(5), BigInteger.valueOf(2) }));
        System.out.println(isIntersect(vertices,low, high));

    }

}
