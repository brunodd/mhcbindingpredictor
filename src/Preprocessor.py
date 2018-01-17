import preprocess.MHCCleaner as MhcCleaner
import preprocess.PeptideCleaner as PeptCleaner
import preprocess.Data as Data

if __name__ == '__main__':
    # Clean MHC data
    MhcCleaner.main()

    # Clean Peptide data
    PeptCleaner.main('data/original/mhc_ligand_full.csv')

    # Combine all data
    Data.main('MhcData.csv', 'AlleleData.csv')