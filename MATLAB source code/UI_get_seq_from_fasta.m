function Seq = UI_get_seq_from_fasta

[file,path] = uigetfile;
if ~isequal(file,0)
    try
        F = fastaread(fullfile(path, file));
        Seq = F(1).Sequence;
    catch
        fig = uifigure;
        MSG = {'The file you chose is either not valid or is not in fasta formant';...
            'You might find it more comfortable to simply paste your sequence into the text box'}
        uialert(fig,MSG,'Invalid File');
        Seq = '';
    end
end
end