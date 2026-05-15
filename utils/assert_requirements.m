function assert_requirements(toolbox_id)
% ASSERT_REQUIREMENTS  Fail fast if MATLAB or add-ons are missing.
%
%   assert_requirements()  — check everything this project needs
%   assert_requirements('Phased_Array_System_Toolbox')
%
%   See also: list_project_requirements

    if nargin < 1
        assert_matlab_version();
        assert_requirements('Phased_Array_System_Toolbox');
        return;
    end

    switch toolbox_id
        case 'Phased_Array_System_Toolbox'
            name = 'Phased Array System Toolbox';
            ok = license('test', 'Phased_Array_System_Toolbox') && ~isempty(ver('phased'));
            if ok
                ok = exist('phased.ULA', 'class') == 8 && exist('physconst', 'file') == 2;
            end
        otherwise
            error('assert_requirements:unknown', 'Unknown toolbox id: %s', toolbox_id);
    end

    if ~ok
        error('pls:requirements:%s', toolbox_id, sprintf([ ...
            '\n*** Additional software required: %s ***\n', ...
            'This project needs the MATLAB add-on above.\n', ...
            'Install: Home tab > Add-Ons > Get Add-Ons > search "%s" > Install\n', ...
            'Then restart MATLAB and run again.\n', ...
            'Dependency map: >> list_project_requirements\n'], ...
            name, name));
    end
end

function assert_matlab_version()
% run_all uses strings(); apply_plot_style uses isgraphics (R2014b+).
    if verLessThan('matlab', '9.1')
        error('pls:requirements:matlab_version', sprintf([ ...
            '\n*** MATLAB R2016b or later required ***\n', ...
            'This project uses string arrays and other features from R2016b+.\n', ...
            'Your release: %s\n'], version));
    end
end
