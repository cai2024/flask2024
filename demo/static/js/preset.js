document.getElementById('WGSBtn').addEventListener('click', function() {initialPreset("WGS",fq_type);});


async function initialPreset(presetType, fq_type) {
    try {
        const response = await fetch('/get_preset', {
            method: 'POST', headers: {'Content-Type': 'application/json'},
            body: JSON.stringify({ presetType: presetType ,mode: fq_type})});
        const data = await response.json();
        
        data.options.forEach(element => {addStep(element)});
        preset_Steps(data.params,data.select)
     	disableAllButtons()
    } catch (error) { console.error('Error initialPreset:', error);}
}


async function preset_Steps(preset_input,preset_select) {
    confirmBtn.disabled = true;  // 禁用确认按钮
    for (let i = 0; i < steps.length; i++) {steps[i].id = i;}
    for (const step of steps) { await preset_OptionsUI(step.id,preset_input,preset_select); } // 初始每一行的选项卡
    $('#selected-steps').sortable('disable'); // 禁用左侧按钮卡排序
    document.getElementById('submitBtn').disabled = false; // 允许提交
}


async function preset_OptionsUI(stepId,preset_input,preset_select) {
    const step = steps.find(step => step.id === stepId);
    const temp_data = await init_type_div(stepId,temp_options) // 获取给定step的输入框+每个选项卡的div内容，获取给定step的推荐选项字典
    document.getElementById('inputArea').insertAdjacentHTML('beforeend', temp_data.div_string); //添加给定step的div
    document.querySelector(`#input-step-${stepId} input`).value = preset_input[stepId]   //读取预设
    step.params = document.querySelector(`#input-step-${stepId} input`).value;      // 初始化step.params
    
    for (let key in preset_select[stepId]) {
        const selectElement = document.querySelector(`#input-step-${stepId} select#select-${stepId}-${key}`)
        const index = Array.from(selectElement.options).findIndex(option => option.value === preset_select[stepId][key]); // 同步推荐
        if (index !== -1) { selectElement.selectedIndex = index; } // 初始化展示的选项
        step.select_option[key] = selectElement.value; }                                    // 初始化step.select_option
    temp_options = await update_my_options(stepId);  //更新提供下一行的总选项
    
    if (step.type.slice(0, 3) === 'map') {
        const index_options =  await index_select(stepId);
        document.getElementById('inputArea').insertAdjacentHTML('beforeend', index_options);
        const selectElement = document.querySelector(`#index-step-${stepId} select`)
        selectElement.onchange = () => { step.select_index=selectElement.value; };
        step.select_index = selectElement.value;                                           // 初始化step.select_index
    }
}
