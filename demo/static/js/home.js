
function toggleImage() {
    var img = document.getElementById('instructionImage');
    img.style.display = (img.style.display === 'none' ? 'block' : 'none');
}




document.addEventListener('keydown', function(event) {
    if (event.ctrlKey && event.key === 'y') {
        event.preventDefault(); // 阻止默认行为，如关闭窗口
        document.addEventListener('keydown', function(e) {
            if (e.key === 'y') {
                document.getElementById('loginPrompt').style.display = 'block';
            }
        }, {once: true});
    }
});



$(document).ready(function() {
    var taskId = $('body').data('task-id');  // 从 HTML body 元素获取 task_id

    $('#folderInput').on('input', function() {
        var query = $(this).val();
        $.ajax({
            url: '/get_folders',
            method: 'GET',
            data: { q: query, task_id: taskId },
            success: function(data) {
                $('#suggestions').empty();
                $.each(data, function(index, folder) {
                    $('#suggestions').append('<div class="folder" data-folder="' + folder + '">' + folder + '</div>');
                });
            }
        });
    });

    $(document).on('click', '.folder', function(e) {
        var folder = $(this).data('folder');
        $('#folderInput').val(folder);
        $('#suggestions').empty();
        $('#folderInput').focus();
        sendFolderSelection(folder);
        e.stopPropagation();
    });

    $(document).on('click', function(e) {
        if (!$(e.target).closest('#suggestions').length && !$(e.target).is('#folderInput')) {
            $('#suggestions').empty();
        }
    });

    $('#folderInput').on('focus', function() {
        var query = $(this).val();
        if (query.trim() !== '') {
            $.ajax({
                url: '/get_folders',
                method: 'GET',
                data: { q: query, task_id: taskId },
                success: function(data) {
                    $('#suggestions').empty();
                    $.each(data, function(index, folder) {
                        $('#suggestions').append('<div class="folder" data-folder="' + folder + '">' + folder + '</div>');
                    });
                }
            });
        }
    });

    function sendFolderSelection(folder) {
        $.ajax({
            url: '/process_folder_selection',
            method: 'POST',
            data: { selected_folder: folder, task_id: taskId },
            success: function(response) {
                console.log("选定的文件夹已成功发送到服务器！");
            },
            error: function(xhr, status, error) {
                console.error("发送文件夹选择时发生错误:", error);
            }
        });
    }
});
