import java.io.*;

import htsjdk.samtools.*;

import java.io.File;
import java.util.*;

public class Parser {
    private String inputReadFile;
    private String inputRefFile;
    private String outputFile;

    public Parser(String inputReadFile, String inputRefFile, String outputFile) {
        this.inputReadFile = inputReadFile;
        this.inputRefFile = inputRefFile;
        this.outputFile = outputFile;

    }

    public void parse() throws IOException {
        // считываем референсный геном
        FileReader reader = new FileReader(this.inputRefFile);
        Scanner scan = new Scanner(reader);
        String ref = scan.nextLine();
        scan.close();
        Map<String, List<Integer>> hashMap = new HashMap<String, List<Integer>>();
        // считываем файл в формате SAM (или BAM)
        File myFile = new File(this.inputReadFile);
        SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(myFile);
        SAMRecordIterator r = sr.iterator();
        int readID = 0;
        while (r.hasNext()) {
            SAMRecord record = r.next();
            //считываем имя рида
            String name = record.getReadName();
            // считываем последовательность нуклеотидов в риде
            String read = record.getReadString();
            // считываем левый индекс выравнивания по референсу
            int leftPosition = record.getAlignmentStart();
            // считываем информацию о мутациях
            String cigar = record.getCigarString();
            readID++;
            this.analyzeCigar(cigar, leftPosition, read, ref, readID, hashMap);

        }
        r.close();
        sr.close();
        this.writer(hashMap);
    }

    private void writer(Map<String, List<Integer>> hashMap) {
        try {
            FileWriter writer = new FileWriter(this.outputFile);
            hashMap.forEach((k, v) -> {
                try {
                    writer.write(k + v.size() + "; " + v + '\n');
                } catch (IOException ex) {
                    System.out.println("Error while writing to output file!");
                }
            });
            writer.close();
        } catch (Exception ex) {
            System.out.println("Error while creating writer");
        }

    }

    private void addToMap(Map<String, List<Integer>> hashMap, String key, int readID) {
        hashMap.computeIfAbsent(key, x -> new ArrayList<>()).add(readID);
    }

    private void analyzeCigar(String cigar, int left, String read, String ref, int readID, Map<String, List<Integer>> hashMap) {
        int refIndex = left - 1;
        int readIndex = 0;
        int shift;
        StringBuilder local = new StringBuilder("");
        for (int i = 0; i < cigar.length(); i++) {
            char c = cigar.charAt(i);
            if (Character.isDigit(c)) {
                local.append(c);
            } else {
                shift = Integer.parseInt(local.toString());
                String key;
                switch (c) {
                    // на этой позиции два нуклеотида
                    case 'M':
                        for (int j = 0; j < shift; j++) {
                            // проверяем, если ли отличающиеся
                            if (read.charAt(readIndex) != ref.charAt(refIndex)) {
                                key = (refIndex) + "; " + ref.charAt(refIndex) + "-->" + read.charAt(readIndex) + "; ";
                                this.addToMap(hashMap, key, readID);
                            }
                            readIndex++;
                            refIndex++;
                        }
                        break;
                    // вставка
                    case 'I':
                        String longInsert = "";
                        for (int j = 0; j < shift; j++) {
                            longInsert += read.charAt(readIndex++);
                        }
                        key = refIndex - 1 + "-" + refIndex + "; " + "_-->" + longInsert + "; ";
                        this.addToMap(hashMap, key, readID);
                        break;
                    // удаление
                    case 'D':
                        for (int j = 0; j < shift; j++) {
                            key = refIndex + "; " + ref.charAt(refIndex++) + "-->_; ";
                            this.addToMap(hashMap, key, readID);
                        }
                        break;
                    // нуклеотиды вырезаны из рида
                    case 'S':
                        readIndex += shift;
                        break;
                    default:
                        System.out.println("Wrong letter in Cigar!");
                }
                local.setLength(0);
            }
        }
    }
}
