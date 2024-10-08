function initializeInputSuggestions(inputSelector, suggestionsId) {
    var inputElement = $(inputSelector);
    if (!inputElement.length) {
        console.error('Specified input element does not exist.');
        return;
    }
    var suggestionsDiv = $('<div>', {
        'id': suggestionsId,
        'class': 'overlay'
    }).insertAfter(inputElement);

    inputElement.on('input focus', function() {
        var folder = $(this).val();
        if (folder.trim() !== '') {
            $.ajax({
                url: '/common/get_folders',
                method: 'GET',
                data: { folder: folder },
                success: function(data) {
                    suggestionsDiv.empty();
                    $.each(data, function(index, folderName) {
                        suggestionsDiv.append($('<div>', {
                            'class': 'folder',
                            'data-folder': folderName,
                            'text': folderName
                        }));
                    });
                }
            });
        } else {
            suggestionsDiv.empty(); // 如果输入为空，清空建议
        }
    });

    $(document).on('click', '.folder', function(e) {
        inputElement.val($(this).data('folder'));
        suggestionsDiv.empty();
        inputElement.focus();
        e.stopPropagation();
    });

    $(document).on('click', function(e) {
        if (!$(e.target).closest('#' + suggestionsDiv.attr('id')).length && !$(e.target).is(inputElement)) {
            suggestionsDiv.empty();
        }
    });
}

function modifyOuterDivClass(inputSelector, newClassName) {
    var outerDiv = $(inputSelector).closest('.form-group.required');
    if (outerDiv.length) {
        outerDiv.addClass(newClassName);
    } else {
        console.error('Outer div not found.');
    }
}
