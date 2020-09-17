import java.io.*;

import htsjdk.samtools.*;
import htsjdk.samtools.*;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import sun.reflect.annotation.ExceptionProxy;

import java.io.File;
import java.io.IOException;
import java.net.InetAddress;
import java.util.*;
import java.util.stream.Collectors;

public class Parser {
    private String inputReadfile;
    private String inputRefFile;
    private String outputFile;

    public Parser(String inputReadfile, String inputRefFile, String outputFile) {
        this.inputReadfile = inputReadfile;
        this.inputRefFile = inputRefFile;
        this.outputFile = outputFile;

    }

    public void parse() throws Exception {
        // считываем референсный геном
        FileReader reader = new FileReader(this.inputRefFile);
        Scanner scan = new Scanner(reader);
        String ref = scan.nextLine();
        Map<String, List<Integer>> hashMap = new HashMap<String, List<Integer>>();
        //System.out.println(this.ref);
        // считываем файл в формате SAM(или BAM)
        File myFile = new File(this.inputReadfile);
        SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(myFile);
        SAMRecordIterator r = sr.iterator();
        int counter = 0;
        int leftPosition;
        String cigar;
        String read;
        while (r.hasNext()) {
            //i++;
            SAMRecord record = r.next();
            //считываем имя рида
            String name = record.getReadName();
            // считываем последовательность нуклеотидов в риде
            read = record.getReadString();
            // считываем левый индекс выравнивания по референсу
            leftPosition = record.getAlignmentStart();
            // считываем информацию о мутациях
            cigar = record.getCigarString();
            counter++;
            this.analyzeCigar(cigar, leftPosition, read, ref, counter, hashMap);
        }
        r.close();
        sr.close();
        //ref = "AAGTCTAGAA";
        //read = "GTCGATAG";
        //cigar = "3M2I3M";
        //leftPosition = 3;
        //this.analyzeCigar(cigar, leftPosition, read, ref, counter, hashMap);
        //System.out.println(hashMap);
        this.writer(hashMap);
    }

    private void writer(Map<String, List<Integer>> hashMap) throws Exception {
        FileWriter writer = new FileWriter(this.outputFile);
        Set<String> keySet = hashMap.keySet();
        for (String key : keySet) {
            List value = hashMap.get(key);
            int length = value.size();
            writer.write(key + length + "; " + value + '\n');
        }
        writer.close();
    }
    private void addToMap(Map<String, List<Integer>> hashMap, String key, int counter){
        if (hashMap.containsKey(key)) {
            hashMap.get(key).add(counter);
        } else {
            List<Integer> list = new ArrayList<Integer>();
            list.add(counter);
            hashMap.put(key, list);
        }
    }

    private void analyzeCigar(String cigar, int left, String read, String ref, int counter, Map<String, List<Integer>> hashMap) {
        int refIndex = left - 1;
        int readIndex = 0;
        int shift;
        String local = "";
        //System.out.println(cigar);
        for (int i = 0; i < cigar.length(); i++) {
            char c = cigar.charAt(i);
            if (Character.isDigit(c)) {
                local = local + c;
            } else {
                shift = Integer.parseInt(local);
                switch (c) {
                    // на этой позиции два нуклеотида
                    case 'M':
                        for (int j = 0; j < shift; j++) {
                            // проверяем, если ли отличающиеся
                            if (read.charAt(j + readIndex) != ref.charAt(j + refIndex)) {
                                String key = refIndex + "; " + ref.charAt(refIndex + j) + "-->" + read.charAt(readIndex + j) + "; ";
                                addToMap(hashMap, key, counter);
                                //System.out.println(refIndex + "; " + ref.charAt(refIndex + j) + "-->" + read.charAt(readIndex + j) + ";");
                            }
                        }
                        refIndex = refIndex + shift;
                        readIndex = readIndex + shift;
                        break;
                    // вставка
                    case 'I':
                        for (int j = readIndex; j < readIndex + shift; j++) {
                            //System.out.println(refIndex - 1 + "-" + refIndex + "; " + "_-->" + read.charAt(j) + ";");
                            String key = refIndex - 1 + "-" + refIndex + "; " + "_-->" + read.charAt(j) + "; ";
                            addToMap(hashMap, key, counter);
                        }
                        readIndex = readIndex + shift;
                        break;
                    // удаление
                    case 'D':
                        for (int j = refIndex; j < refIndex + shift; j++) {
                            //System.out.println(j + "; " + ref.charAt(j) + "-->_;");
                            String key = j + "; " + ref.charAt(j) + "-->_; ";
                            addToMap(hashMap, key, counter);
                        }
                        refIndex = refIndex + shift;
                        break;
                    // нуклеотиды вырезаны из рида
                    case 'S':
                        readIndex = readIndex + shift;
                        break;
                    default:
                        System.out.println("Wrong letter in Cigar!");
                }
                local = "";
            }
        }
    }
}