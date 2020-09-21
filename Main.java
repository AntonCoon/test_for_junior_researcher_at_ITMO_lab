import java.io.IOException;

public class Main {
    public static void main(String args[]) {
        if (args.length != 3){
            System.out.println("You need three arguments: .sam (or bam) .txt(reference) .txt(output) ");
            System.exit(1);
        }
        Parser parser = new Parser(args[0], args[1], args[2]);
        try {
            parser.parse();
        } catch (IOException e) {
            System.out.println("IOException. You need three arguments: .smm (or bmm) .txt(reference) .txt(output) ");
        }
    }
}
