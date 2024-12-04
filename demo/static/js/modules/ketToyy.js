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

